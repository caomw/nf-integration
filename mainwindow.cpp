#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <QFileDialog>

#include <iostream>

using namespace std;

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow),
    m_imgs(),
    m_mesh(),
    m_mesh_filename(),
    m_target_normal(),
    m_viewer(new CTriMeshViewer(&m_mesh,parent)) {

    ui->setupUi(this);

    connect(this,SIGNAL(cameraChanged(CCamera<float>)),m_viewer,SLOT(updateCam(CCamera<float>)));
    connect(this,SIGNAL(viewpointChanged(CViewPoint<float>)),m_viewer,SLOT(updateView(CViewPoint<float>)));

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
    if(!m_mesh.has_vertex_colors() || !m_mesh.has_face_colors()) {

        m_mesh.request_vertex_colors();
        m_mesh.request_face_colors();

    }

    // set vertex/face colors to default
    TriangleMesh::VertexIter v_it;
    for (v_it = m_mesh.vertices_begin(); v_it != m_mesh.vertices_end(); ++v_it)
        m_mesh.set_color(v_it,OpenMesh::Vec3uc(200,200,200));
    TriangleMesh::FaceIter f_it;
    for (f_it = m_mesh.faces_begin(); f_it != m_mesh.faces_end(); ++f_it)
        m_mesh.set_color(f_it,OpenMesh::Vec3uc(200,200,200));

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

    // set maximum value of spin box and trigger redraw
    ui->imgSpinBox->setMaximum(m_imgs.size()-1);
    ui->imgSpinBox->setValue(m_imgs.size()-1);
    emit viewpointChanged(m_imgs.back().GetViewpoint());
    emit cameraChanged(m_imgs.back().GetCam());
    m_viewer->update();

}

void MainWindow::on_actionLoad_Mesh_triggered() {

    m_mesh_filename = QFileDialog::getOpenFileName(this, tr("Open start mesh..."),
                                                   ".",
                                                   tr("(*.ply);;(*.stl);;(*.off);;(*.obj)"));

    // everything else is the same as in reload
    on_actionReload_triggered();

    // first load
    m_viewer->updateBoundingBox();

}

void MainWindow::on_imgSpinBox_valueChanged(int arg1) {

    if(arg1>=0 && arg1<m_imgs.size() && m_imgs.size()>0) {

        emit cameraChanged(m_imgs.at(arg1).GetCam());
        emit viewpointChanged(m_imgs.at(arg1).GetViewpoint());

    }

}
