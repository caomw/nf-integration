#include <iostream>
#include <QFileDialog>

#include "viewer.h"

#define NEAR_PLANE_TOL 0.9999

using namespace std;

CViewer::CViewer(QWidget* parent):
    QGLWidget(parent),
    m_cam(),
    m_viewpoint(),
    m_znear(0.1),
    m_zfar(50.0),
    m_last_point(),
    m_center(),
    m_bbox(),
    m_show_color(true) {

    // set up window size according to resolution of camera
    CVector<size_t,2> sizes = m_cam.GetSize();
    this->setFixedSize(sizes.Get(0),sizes.Get(1));

    // set cursor type and window title
    setCursor(Qt::PointingHandCursor);
    setWindowTitle("OpenGL Viewer");

}

CViewer::CViewer(const CCamera<float>& cam, QWidget* parent):
    QGLWidget(parent),
    m_cam(cam),
    m_viewpoint(),
    m_znear(0.1),
    m_zfar(50.0),
    m_last_point(),
    m_center(),
    m_bbox(),
    m_show_color(true) {

    // set up window size according to resolution of camera
    CVector<size_t,2> sizes = m_cam.GetSize();
    this->setFixedSize(sizes.Get(0),sizes.Get(1));

    // set cursor type and window title
    setCursor(Qt::PointingHandCursor);
    setWindowTitle("OpenGL Viewer");

}

void CViewer::updateCam(const CCamera<float>& cam) {

    m_cam = cam;

    CVector<size_t,2> sizes = m_cam.GetSize();
    this->setFixedSize(sizes.Get(0),sizes.Get(1));

    glViewport(0,0,sizes.Get(0),sizes.Get(1));

    loadProjectionMatrix();

    this->updateGL();

}

void CViewer::loadProjectionMatrix() {

    // get camera intrinsics
    matf K = m_cam.GetOpenGLProjectionMatrix(m_znear,m_zfar);

    float mp[16];

    // dense matrix stores in row-major order now!!!
    mp[0] = K.Get(0,0);
    mp[1] = K.Get(1,0);
    mp[2] = K.Get(2,0);
    mp[3] = K.Get(3,0);
    mp[4] = K.Get(0,1);
    mp[5] = K.Get(1,1);
    mp[6] = K.Get(3,1);
    mp[7] = K.Get(3,1);
    mp[8] = K.Get(0,2);
    mp[9] = K.Get(1,2);
    mp[10] = K.Get(2,2);
    mp[11] = K.Get(3,2);
    mp[12] = K.Get(0,3);
    mp[13] = K.Get(1,3);
    mp[14] = K.Get(2,3);
    mp[15] = K.Get(3,3);

    glMatrixMode(GL_PROJECTION);
    glLoadMatrixf(&mp[0]);

}

void CViewer::loadView(const CRigidMotion<float,3>& viewpoint) {

    matf F = matf(viewpoint);
    float mv[16];

    // dense matrix stores in row-major order!!!
    mv[0] = F.Get(0,0);
    mv[1] = -F.Get(1,0);
    mv[2] = -F.Get(2,0);
    mv[3] = F.Get(3,0);
    mv[4] = F.Get(0,1);
    mv[5] = -F.Get(1,1);
    mv[6] = -F.Get(2,1);
    mv[7] = F.Get(3,1);
    mv[8] = F.Get(0,2);
    mv[9] = -F.Get(1,2);
    mv[10] = -F.Get(2,2);
    mv[11] = F.Get(3,2);
    mv[12] = F.Get(0,3);
    mv[13] = -F.Get(1,3);
    mv[14] = -F.Get(2,3);
    mv[15] = F.Get(3,3);

    glMatrixMode(GL_MODELVIEW);
    glLoadMatrixf(&mv[0]);

}

void CViewer::initializeGL(){

    glClearColor(0.1,0.1,0.1,1.0);

    // this is in normalized device coordinates
    glClearDepth(1.0);

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_COLOR_MATERIAL);
    glDepthFunc(GL_LEQUAL);

    // set viewport size != window size
    QSize ws = this->size();
    glViewport(0,0,ws.width(),ws.height());

    // load intrinsics
    loadProjectionMatrix();

    // load extrinsics
    loadView(m_viewpoint.GetTransformation());

    // set light positions
    GLfloat pos1[] = { 0.1,  0.1, -0.02, 0.0};
    GLfloat pos2[] = {-0.1,  0.1, -0.02, 0.0};
    GLfloat pos3[] = { 0.0,  0.0,  0.1,  0.0};
    GLfloat col1[] = { 0.7,  0.7,  0.8,  1.0};
    GLfloat col2[] = { 0.8,  0.7,  0.7,  1.0};
    GLfloat col3[] = { 0.5,  0.5,  0.5,  1.0};

    glEnable(GL_LIGHT0);
    glLightfv(GL_LIGHT0,GL_POSITION, pos1);
    glLightfv(GL_LIGHT0,GL_DIFFUSE,  col1);
    glLightfv(GL_LIGHT0,GL_SPECULAR, col1);

    glEnable(GL_LIGHT1);
    glLightfv(GL_LIGHT1,GL_POSITION, pos2);
    glLightfv(GL_LIGHT1,GL_DIFFUSE,  col2);
    glLightfv(GL_LIGHT1,GL_SPECULAR, col2);

    glEnable(GL_LIGHT2);
    glLightfv(GL_LIGHT2,GL_POSITION, pos3);
    glLightfv(GL_LIGHT2,GL_DIFFUSE,  col3);
    glLightfv(GL_LIGHT2,GL_SPECULAR, col3);

    glEnable(GL_LIGHTING);
    glShadeModel(GL_SMOOTH);

    //int zbuffer_depth;
    //glGetIntegerv(GL_DEPTH_BITS, &zbuffer_depth);
    //cout << "Buffer depth: " << zbuffer_depth << endl;

}

void CViewer::paintGL(){

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // paint some dummy object
    glBegin(GL_TRIANGLES);

    if(m_show_color)
        glColor3f(0.0f,0.0f,1.0f);
    else
        glColor3f(0.5f,0.5f,0.5f);

    glVertex3f(0.0f,0.0f,2.0f);

    if(m_show_color)
        glColor3f(0.0f,1.0f,0.0f);

    glVertex3f(1.0f,0.0f,2.0f);

    if(m_show_color)
        glColor3f(1.0f,0.0f,0.0f);

    glVertex3f(0.0f,1.0f,2.0f);

    glEnd();

    // draw some points
    glPointSize(5.0f);
    glBegin(GL_POINTS);
    glColor3f(1.0f,1.0f,1.0f);
    glVertex3f(0.0f,2.0f,2.0f);
    glColor3f(1.0f,1.0f,1.0f);
    glVertex3f(2.0f,2.0f,2.0f);
    glColor3f(1.0f,1.0f,1.0f);
    glVertex3f(2.0f,0.0f,2.0f);
    glEnd();

}

void CViewer::wheelEvent(QWheelEvent* event) {

    float radius = 1;
    float d = -float(event->delta())/120.0*0.2*radius;

    CVector<float,3> dt = { 0, 0, 0.01f * (m_zfar - m_znear) * d };
    m_viewpoint.DifferentialTranslate(dt);

    this->loadView(m_viewpoint.GetTransformation());
    this->updateGL();

    event->accept();

}

void CViewer::mousePressEvent(QMouseEvent* event) {

    m_last_point = event->pos();
    event->accept();

}

void CViewer::mouseMoveEvent(QMouseEvent* event) {

    QPoint newpoint = event->pos();

    float dx = newpoint.x() - m_last_point.x();
    float dy = newpoint.y() - m_last_point.y();

    if (event->buttons()==Qt::MidButton) {

        // need depth here to backproject image plane translations (take mean of clipping depths)
        CVector<float,3> dt = { -0.0005f*(m_zfar-m_znear)*dx, -0.0005f*(m_zfar-m_znear)*dy, 0.0f };
        m_viewpoint.DifferentialTranslate(dt);

    }
    else if (event->buttons() == Qt::LeftButton)
    {

        CVector<float,3> axis = { -0.05f*dy, 0.05f*dx, 0 };

        // transform rotation axis into world coordinates
        const CRigidMotion<float,3> Finv = m_viewpoint.GetInverseTransformation();
        axis = Finv.DifferentialTransform(axis);

        m_viewpoint.Orbit(m_center,axis);

    }
    else if (event->buttons() == Qt::RightButton)
    {

        CVector<float,3> dt;
        if(dx<0)
            dt(2) = 0.0005*(m_zfar-m_znear)*fabs(dx);
        else
            dt(2) = -0.0005*(m_zfar-m_znear)*fabs(dx);

        m_viewpoint.DifferentialTranslate(dt);

    }

    this->loadView(m_viewpoint.GetTransformation());
    this->updateGL();

    // remember this point
    m_last_point = newpoint;

    event->accept();

}

void CViewer::keyPressEvent(QKeyEvent* event) {

    // reset to cam to world coordinates
    if(event->key() == Qt::Key_Z) {

        CRigidMotion<float,3> id;
        this->loadView(id);

        this->updateGL();

    }

    if(event->key() == Qt::Key_Minus) {

        CVector<float,3> dt;
        dt(2) = -0.05*(m_zfar-m_znear);
        m_viewpoint.DifferentialTranslate(dt);
        this->loadView(m_viewpoint.GetTransformation());
        this->updateClipDepth(m_viewpoint.GetTransformation(),NEAR_PLANE_TOL);
        this->updateGL();

    }

    if(event->key() == Qt::Key_Plus) {

        CVector<float,3> dt;
        dt(2) = 0.05*(m_zfar-m_znear);
        m_viewpoint.DifferentialTranslate(dt);
        this->loadView(m_viewpoint.GetTransformation());
        this->updateClipDepth(m_viewpoint.GetTransformation(),NEAR_PLANE_TOL);
        this->updateGL();

    }

    // turn color rendering on/off
    if(event->key() == Qt::Key_C) {

        m_show_color = !m_show_color;
        this->updateGL();

    }

    // save screenshot
    if(event->key() == Qt::Key_S) {

        QString filename = QFileDialog::getSaveFileName(this, tr("Export screenshot..."),
                                                        ".",
                                                        tr("(*.png)"));

        QImage img = this->grabFrameBuffer();
        img.save(filename);

    }


}

CDenseArray<float> CViewer::getDepthMap(const CRigidMotion<float,3>& viewpoint) {

    // load view and paint
    this->loadView(viewpoint);
    this->updateClipDepth(viewpoint,NEAR_PLANE_TOL);
    this->paintGL();

    // allocate result
    QSize ws = this->size();
    matf z(ws.height(),ws.width());

    cout << viewpoint << endl;
    //cout << "Bounding box: " << endl;
    //cout << m_bbox.Lower() << " " << m_bbox.Upper() << endl;
    //cout << "Clip depths: " << endl;
    //cout << m_znear << " " << m_zfar << endl;

    // we need buffer in order to flip the array
    float* buffer = new float[z.NElems()];

    // read depths
    glReadPixels(0, 0,z.NCols(),z.NRows(),GL_DEPTH_COMPONENT,GL_FLOAT,buffer);

    // because we have updated the depths at the beginning, we don't need to grab them from the graphics card
    float zfmzn = m_zfar - m_znear;
    float zfpzn = m_zfar + m_znear;
    float zftzn = 2.0*m_zfar*m_znear;

    for(size_t i=0; i<z.NRows(); i++) {

        for(size_t j=0; j<z.NCols(); j++) {

            // normalized depth
            float zn = 2.0*buffer[(z.NRows()-i-1)*z.NCols() + j] - 1.0;

            z(i,j) = zftzn/(zfpzn - zn*zfmzn);

        }


    }

    delete [] buffer;

    // restore view and repaint
    this->loadView(m_viewpoint.GetTransformation());
    this->updateClipDepth(m_viewpoint.GetTransformation(),NEAR_PLANE_TOL);
    this->paintGL();

    return z;

}


CTriMeshViewer::CTriMeshViewer(const CCamera<float>& cam, const CTriangleMesh* mesh, QWidget* parent):
    CViewer(cam,parent),
    m_mesh(mesh) {

    if(m_mesh!=nullptr)
        this->updateBoundingBox();

}

void CTriMeshViewer::updateBoundingBox() {

    if(m_mesh==nullptr)
        return;

    // recompute bounding box
    m_bbox = m_mesh->BoundingBox();

    // set center of camera rotation
    m_center = m_bbox.Barycenter();

    // update clip depths
    this->updateClipDepth(m_viewpoint.GetTransformation(),NEAR_PLANE_TOL);

}

void CTriMeshViewer::updateView(const CViewPoint<float>& viewpoint) {

    // set the member variable
    m_viewpoint = viewpoint;

    // send the view to the graphics card
    this->loadView(viewpoint.GetTransformation());

    // update clip depths
    this->updateClipDepth(viewpoint.GetTransformation());

    // render
    this->updateGL();

}

void CTriMeshViewer::updateClipDepth(const CRigidMotion<float,3>& F, float tolerance) {

    if(m_mesh==nullptr)
        return;

    vector<CVector<float,3> > corners = m_bbox.Corners();

    float minz = std::numeric_limits<float>::max();
    float maxz = -std::numeric_limits<float>::max();

    for(u_int i=0; i<corners.size(); i++) {

        CVector<float,3> lc = F.Transform(corners.at(i));

        if(lc.Get(2)<minz)
            minz = lc.Get(2);

        if(lc.Get(2)>maxz)
            maxz = lc.Get(2);

    }

    float ntol, ftol;

    if(tolerance==0)
        ntol = 1.0;
    else
        ntol = tolerance;

    ftol = 1.0/ntol;

    m_znear = ntol*minz;

    if(m_znear<=0)
        m_znear = 1e-3;

    float newfar = ftol*maxz;
    if(newfar>m_znear)
        m_zfar = newfar;

    // send new projection matrix to graphics card
    this->loadProjectionMatrix();

}

void CTriMeshViewer::mouseReleaseEvent(QMouseEvent *event) {

    // update clip depths only after change is done
    this->updateClipDepth(m_viewpoint.GetTransformation(),NEAR_PLANE_TOL);
    this->updateGL();
    event->accept();

}

void CTriMeshViewer::paintGL() {

    if(m_mesh==nullptr)
        return;

    // set projection and modelview matrices from m_view
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // load mesh
    TriangleMesh::ConstFaceIter f_it;
    TriangleMesh::ConstFaceVertexIter fv_it;

    glBegin(GL_TRIANGLES);

    for (f_it=m_mesh->faces_begin(); f_it!=m_mesh->faces_end(); ++f_it) {

        // load normal, FIXME: are the normals computed in the wrong way?
        // or does OpenGL flip them???
        TriangleMesh::Normal temp = -m_mesh->normal(f_it);
        glNormal3fv(&temp[0]);

        if(m_show_color) {

            OpenMesh::Vec3uc color = m_mesh->color(f_it);
            glColor3ub(color[0],color[1],color[2]);

        }
        else
            glColor3ub(200,200,200);

        for (fv_it = m_mesh->cfv_iter(f_it.handle()); fv_it; ++fv_it) {

            // load vertex
            glVertex3fv(&m_mesh->point(fv_it)[0]);

        }

    }

    glEnd();

}

