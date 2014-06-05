#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "flowviz.h"

#include <QFileDialog>

#include <iostream>
#import <chrono>

using namespace std;

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow),
    m_imgs(),
    m_mesh(),
    m_mesh_filename(),
    m_viewer(new CTriMeshViewer(&m_mesh,parent)) {

    ui->setupUi(this);

    connect(this,SIGNAL(cameraChanged(CCamera<float>)),m_viewer,SLOT(updateCam(CCamera<float>)));
    connect(this,SIGNAL(viewpointChanged(CViewPoint<float>)),m_viewer,SLOT(updateView(CViewPoint<float>)));
    connect(this->ui->imgSpinBox,SIGNAL(valueChanged(int)),this,SLOT(showNormalFieldImage(int)));

    OpenMesh::FPropHandleT<vec3f> nt;
    m_mesh.add_property(nt,"nt");

}

MainWindow::~MainWindow()
{

    delete m_viewer;
    delete ui;

}

void MainWindow::on_actionExit_triggered() {

    QApplication::exit();

}

void MainWindow::on_action3D_Scene_triggered() {

    m_viewer->show();

}

void MainWindow::on_actionReload_triggered() {

    if(m_mesh_filename.isEmpty())
        return;

    if(!OpenMesh::IO::read_mesh(m_mesh,m_mesh_filename.toStdString().c_str())) {

        ui->statusBar->showMessage("Could not load mesh...");
        return;

    }

    // make sure we have normals
    if(!m_mesh.has_vertex_normals() || !m_mesh.has_face_normals()) {

         m_mesh.request_face_normals();
         m_mesh.request_vertex_normals();

    }

    m_mesh.update_face_normals();
    m_mesh.update_vertex_normals();
    m_mesh.update_normals();

    // attach color
    if(!m_mesh.has_face_colors())
        m_mesh.request_face_colors();

    // set vertex/face colors to default
    TriangleMesh::FaceIter f_it;
    for (f_it = m_mesh.faces_begin(); f_it != m_mesh.faces_end(); ++f_it)
        m_mesh.set_color(f_it,OpenMesh::Vec3uc(200,200,200));

    m_viewer->updateBoundingBox();
    m_viewer->updateGL();

}

void MainWindow::on_actionOpen_View_triggered() {

    bool wasempty = m_imgs.empty();

    QStringList filenames = QFileDialog::getOpenFileNames(this,tr("Open images..."),".",tr("*.exr"));

    if(filenames.size()==0)
        return;

    for(size_t i=0; i<filenames.size(); i++) {

        try {

            CNormalField<float> nf;
            nf.ReadFromFile(filenames.at(i).toStdString().c_str());
            m_imgs.push_back(nf);

        }
        catch(const std::runtime_error& e) {
            ui->statusBar->showMessage(e.what());
        }

    }

    // if these are the first images, enable navigation elements in gui
    if(wasempty)
        ui->imgSpinBox->setEnabled(true);

    if(m_imgs.size()==1) {

        this->on_imgSpinBox_valueChanged(0);
        this->showNormalFieldImage(0);

    }
    else {

        ui->imgSpinBox->setMaximum(m_imgs.size()-1);
        ui->imgSpinBox->setValue(m_imgs.size()-1);

    }
    //this->on_imgSpinBox_valueChanged(m_imgs.size()-1);
    //emit viewpointChanged(m_imgs.back().GetViewpoint());
    //emit cameraChanged(m_imgs.back().GetCam());

    m_viewer->update();

}

void MainWindow::on_actionLoad_Mesh_triggered() {

    m_mesh_filename = QFileDialog::getOpenFileName(this, tr("Open start mesh..."),
                                                   ".",
                                                   tr("(*.ply);;(*.stl);;(*.off);;(*.obj)"));

    // everything else is the same as in reload
    on_actionReload_triggered();

}

void MainWindow::on_imgSpinBox_valueChanged(int arg1) {

    if(arg1>=0 && arg1<m_imgs.size() && m_imgs.size()>0) {

        emit cameraChanged(m_imgs.at(arg1).GetCam());
        emit viewpointChanged(m_imgs.at(arg1).GetViewpoint());
        cout << "bla" << endl;
    }


}

void MainWindow::showNormalFieldImage(int no) {

    CFlowVisualization<float> vis;
    QImage img = vis.CalcDirectionEncoding(m_imgs.at(no));
    QSize slabel = ui->imgLabel->size();
    ui->imgLabel->setPixmap(QPixmap::fromImage(img.scaled(slabel)));

}

void MainWindow::on_stepButton_clicked() {

    chrono::high_resolution_clock::time_point t0, t1;
    chrono::duration<double> time_span;

    // compute residual
    CDenseVector<double> r;
    double total_error = computeResidual(r);
    stringstream ss;
    ss << "Total error: " << total_error;
    ui->statusBar->showMessage(ss.str().c_str());

    // assemble discrete gradient operator
    t0 = chrono::high_resolution_clock::now();
    CCSCMatrix<double,int> G = m_mesh.ComputeGradientOperator();
    t1 = chrono::high_resolution_clock::now();
    time_span = chrono::duration_cast<chrono::duration<double>>(t1 - t0);
    cout << "Time assembly: " << time_span.count() << " seconds." << endl;
    G.SaveToFile("/home/jbalzer/G.txt");
    // compute right-hand side, but make sure that b accounts for normalization condition
    CDenseVector<double> b(G.NCols()+1);
    G.Multiply(b,r);

    // square matrix, resize, and add normalization
    t0 = chrono::high_resolution_clock::now();
    CCSCMatrix<double,int> GtG = CCSCMatrix<double,int>::Square(G);
    t1 = chrono::high_resolution_clock::now();
    time_span = chrono::duration_cast<chrono::duration<double>>(t1 - t0);
    cout << "Time squaring: " << time_span.count() << " seconds."  << endl;


    GtG.Resize(m_mesh.n_vertices()+1,m_mesh.n_vertices()+1);
    size_t nnz_old = GtG.NNz();

    vector<int>* rowptr = GtG.GetRowIndices().get();
    vector<int>* colptr = GtG.GetColumnPointer().get();
    vector<double>* valptr = GtG.GetValues().get();

    for(size_t k=0; k<m_mesh.n_vertices(); k++) {
        rowptr->push_back(k);
        valptr->push_back(1.0);
    }
    rowptr->push_back(m_mesh.n_vertices());
    valptr->push_back(-1.0);
    colptr->at(m_mesh.n_vertices()+1) = nnz_old + m_mesh.n_vertices()+1;

    // FIXME: Reuse factor!!!
    CCholeskySolver<double> solver(GtG);
    CDenseArray<double> vn = solver.Solve(GtG,b);

    // deform the mesh mesh, use one vertex less
    for(size_t i=0; i<vn.NElems()-1; i++) {
      OpenMesh::VertexHandle vh(i);
      m_mesh.point(vh) = m_mesh.point(vh) + m_mesh.normal(vh)*vn.Get(i,0);

    }

    // update normals
    m_mesh.update_normals();


    m_viewer->updateGL();

}

double MainWindow::computeResidual(CDenseVector<double>& r) {

    // FIXME: compute visible views

    if(r.NElems()!=3*m_mesh.n_faces())
        r = vecf(3*m_mesh.n_faces());

    double total_error = 0;
    TriangleMesh::FaceIter f_it;
    size_t row = 0;

    for(f_it=m_mesh.faces_begin(); f_it!=m_mesh.faces_end(); ++f_it) {

        // barycenter of face normal
        const vec3f x = m_mesh.Barycenter(f_it);
        const vec3f n = m_mesh.Normal(f_it);

        // compute normal estimate
        vec3f nt;
        switch(ui->intSourceComboBox->currentIndex()) {

            case 0:
                nt = m_imgs.back().Get(x);
                break;
            case 1:
                nt = m_imgs.back().GetDeflectometricNormal(x);
                break;

        }

        // normalize
        if(nt.Norm2()>0)
            nt.Normalize();

        OpenMesh::FPropHandleT<vec3f> ntp;

        if(m_mesh.get_property_handle(ntp,"nt"))
            m_mesh.property(ntp,f_it) = nt;

        // error
        vec3f dn = n - nt;

        // norm of error, between 0,2
        float rn;
        if(nt.Norm2()>0)
            rn = dn.Norm2()/2.0;
        else
            rn = 0;

        total_error += rn;

        // set color
        QColor val = QColor::fromHsvF((1.0-rn)*0.6667,1.0,1.0);
        m_mesh.set_color(f_it,OpenMesh::Vec3uc(val.red(),val.green(),val.blue()));

        // quadrature weight
        float A =  m_mesh.FaceArea(f_it);
        A = sqrt(A);

        for(uint i=0; i<3; i++)
            r(row+i) = dn.Get(i)*A;

        // increment row pointer
        row += 3;

    }

    return total_error;

}

void MainWindow::on_actionSave_Mesh_triggered() {

    QString type; // = 0

    QString filename = QFileDialog::getSaveFileName(this, tr("Export mesh..."),
                                                    ".",
                                                    tr("*.ply;;*.vtk;;*.stl;;*.off;;*.obj"),
                                                    &type);

    if(!filename.endsWith(".off") && !filename.endsWith(".stl") && !filename.endsWith(".obj") && !filename.endsWith(".ply") && !filename.endsWith(".vtk"))
        filename += type.remove(0,1);

    try {

        m_mesh.SaveToFile(filename.toStdString().c_str());

    }
    catch(const std::runtime_error& e) {

        ui->statusBar->showMessage(e.what());

    }

}

void MainWindow::on_actionError_triggered() {

    CDenseVector<double> r;
    double total_error = computeResidual(r);
    stringstream ss;
    ss << "Total error: " << total_error;
    ui->statusBar->showMessage(ss.str().c_str());
    m_viewer->updateGL();

}

void MainWindow::on_actionRefine_triggered() {

    m_mesh.UniformMeshRefinement(1);

    // set default vertex/face color
    TriangleMesh::FaceIter f_it;
    for (f_it = m_mesh.faces_begin(); f_it != m_mesh.faces_end(); ++f_it)
        m_mesh.set_color(f_it,OpenMesh::Vec3uc(200,200,200));


    m_mesh.update_normals();
    m_viewer->updateGL();

}

void MainWindow::on_actionSave_Image_triggered() {

    if(m_imgs.size()==0)
        return;

    CFlowVisualization<float> vis;
    QImage img = vis.CalcDirectionEncoding(m_imgs.at(ui->imgSpinBox->value()));

    QString type;
    QString filename = QFileDialog::getSaveFileName(this, tr("Export mesh..."),
                                                    ".",
                                                    tr("*.png;;*.jpg"),
                                                    &type);

    if(!filename.endsWith(".png") && !filename.endsWith(".jpg"))
        filename += type.remove(0,1);

    img.save(filename);

}
