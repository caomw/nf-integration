#include <iostream>

#include "viewer.h"

#define NEAR_PLANE_TOL 0.8
#define FAR_PLANE_TOL 1.6

using namespace std;

CViewer::CViewer(const CPinholeCam<float>& cam, QWidget* parent):
    QGLWidget(parent),
    m_cam(cam),
    m_viewpoint(),
    m_znear(0.1),
    m_zfar(100.0),
    m_last_point(),
    m_center(),
    //m_bbox(),
    m_show_color(true) {

    // set up window size according to resolution of camera
    CVector<size_t,2> sizes = m_cam.GetSize();
    setFixedSize(sizes.Get(0),sizes.Get(1));

    // set cursor type and window title
    setCursor(Qt::PointingHandCursor);
    setWindowTitle("OpenGL Viewer");

}

void CViewer::loadProjectionMatrix() {

    // get camera intrinsics
    matf K = m_cam.GetOpenGLProjectionMatrix(m_znear,m_zfar);
    glLoadMatrixf(K.Data().get());

}

void CViewer::loadView(const CRigidMotion<float,3>& viewpoint) {

    matf F = matf(viewpoint);

    float mv[16];
    mv[0] = F.Data().get()[0];
    mv[1] = -F.Data().get()[1];
    mv[2] = -F.Data().get()[2];
    mv[3] = F.Data().get()[3];
    mv[4] = F.Data().get()[4];
    mv[5] = -F.Data().get()[5];
    mv[6] = -F.Data().get()[6];
    mv[7] = F.Data().get()[7];
    mv[8] = F.Data().get()[8];
    mv[9] = -F.Data().get()[9];
    mv[10] = -F.Data().get()[10];
    mv[11] = F.Data().get()[11];
    mv[12] = F.Data().get()[12];
    mv[13] = -F.Data().get()[13];
    mv[14] = -F.Data().get()[14];
    mv[15] = F.Data().get()[15];

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

    loadView(m_viewpoint.GetTransformation());
    updateGL();

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
        CVector<float,3> dt = { -0.0005*(m_zfar-m_znear)*dx, -0.0005*(m_zfar-m_znear)*dy, 0 };
        m_viewpoint.DifferentialTranslate(dt);

    }
    else if (event->buttons() == Qt::LeftButton)
    {

        CVector<float,3> axis = { -0.05*dy, 0.05*dx, 0 };

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

    loadView(m_viewpoint.GetTransformation());
    updateGL();

    // remember this point
    m_last_point = newpoint;

    event->accept();

}

void CViewer::keyPressEvent(QKeyEvent* event) {

    // reset to cam to world coordinates
    if(event->key() == Qt::Key_Z) {

        CRigidMotion<float,3> id;
        loadView(id);

        updateGL();

    }

    if(event->key() == Qt::Key_B) {

        updateBoundingBox();

        //m_center = m_bbox.Barycenter();

        updateGL();

    }

    // turn color rendering on/off
    if(event->key() == Qt::Key_C) {

        m_show_color = !m_show_color;

        updateGL();

    }

}

CDenseArray<float> CViewer::getDepthMap(const CRigidMotion<float,3>& viewpoint) {

    // load view and paint
    loadView(viewpoint);
    updateClipDepth(viewpoint,1.0);
    paintGL();

    // allocate result
    QSize ws = this->size();
    matf z(ws.height(),ws.width());

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
    loadView(m_viewpoint.GetTransformation());
    updateClipDepth(m_viewpoint.GetTransformation(),NEAR_PLANE_TOL);
    paintGL();

    return z;

}


/*CTriMeshViewer::CTriMeshViewer(const CView<float>& view, const CTriangleMesh* mesh, QWidget* parent):
    CViewer(view,parent),
    m_mesh(mesh) {

    if(m_mesh!=nullptr)
        updateBoundingBox();

}

void CTriMeshViewer::setMesh(const CTriangleMesh *mesh) {

    m_mesh = mesh;

    // update bounding box
    updateBoundingBox();

}

void CTriMeshViewer::updateBoundingBox() {

    if(m_mesh==nullptr)
        return;

    // recompute bounding box
    m_bbox = m_mesh->BoundingBox();

    // update clip depths
    updateClipDepth(m_view,NEAR_PLANE_TOL);

}

CDenseArray<int> CTriMeshViewer::getFaceMap(const CView<float>& view) {

    // allocate result
    QSize ws = this->size();
    CDenseArray<int> fm(ws.height(),ws.width());

    if(m_mesh==nullptr)
        return fm;

    // turn lighting off and set clear color to black
    glDisable(GL_LIGHTING);
    glClearColor(0.0,0.0,0.0,1.0);

    // load view
    loadView(view);

    // clear buffers
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // load mesh
    TriangleMesh::ConstFaceIter f_it;
    TriangleMesh::ConstFaceVertexIter fv_it;

    //unsigned int counter = 1;

    // draw it
    glBegin(GL_TRIANGLES);

    for (f_it=m_mesh->faces_begin(); f_it!=m_mesh->faces_end(); ++f_it) {

        unsigned int id = f_it.handle().idx() + 1;

        // take the last three bytes of the counter as the color
        uchar* bytes = reinterpret_cast<uchar*>(&id);

        // load color
        glColor3ub(bytes[0],bytes[1],bytes[2]);

        // draw vertices
        fv_it = m_mesh->cfv_iter(f_it.handle());
        glVertex3fv(&m_mesh->point(fv_it)[0]);
        ++fv_it;
        glVertex3fv( &m_mesh->point(fv_it)[0] );
        ++fv_it;
        glVertex3fv( &m_mesh->point(fv_it)[0] );

    }

    glEnd();

    glFlush();

    // need buffer in order to flip the image
    uchar* buffer = new uchar[3*fm.NElems()];

    glReadPixels(0, 0,fm.NCols(),fm.NRows(),GL_RGB,GL_UNSIGNED_BYTE,buffer);

    for(size_t i=0; i<fm.NRows(); i++) {

        for(size_t j=0; j<fm.NCols(); j++) {

            // normalized depth
            uchar* pbuffer = buffer + 3*((fm.NRows()-i-1)*fm.NCols() + j);

            // concatenate the three bytes to an integer
            unsigned int index = 0;
            uchar* bytes = reinterpret_cast<uchar*>(&index);
            bytes[0] = pbuffer[0];
            bytes[1] = pbuffer[1];
            bytes[2] = pbuffer[2];

            // transform to camera coordinates
            if(index>0)
                fm(i,j) = index - 1;
            else
                fm(i,j) = -1;

        }

    }

    // clean up
    delete [] buffer;

    // restore view and graphics display
    glEnable(GL_LIGHTING);
    glClearColor(0.1,0.1,0.1,1.0);
    loadView(m_view);
    paintGL();

    return fm;

}

void CTriMeshViewer::updateView(const R4R::CView<float>& view) {

    // set the member variable
    m_view = view;

    // send the view to the graphics card
    loadView(view);

    // update clip depths
    updateClipDepth(view);

    // render
    updateGL();

}

void CTriMeshViewer::updateClipDepth(const CView<float>& view, float tolerance) {

    if(m_mesh==nullptr)
        return;

    vector<CVector<float,3> > corners = m_bbox.Corners();

    float minz = std::numeric_limits<float>::max();
    float maxz = -std::numeric_limits<float>::max();

    CRigidMotion<float,3> F = view.GetTransformation();

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
    loadProjectionMatrix();

}

void CTriMeshViewer::mouseReleaseEvent(QMouseEvent *event) {

    // update clip depths only after change is done
    this->updateClipDepth(m_view,NEAR_PLANE_TOL);

    updateGL();

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

        // load normal
        glNormal3fv(&m_mesh->normal(f_it)[0]);

        for (fv_it = m_mesh->cfv_iter(f_it.handle()); fv_it; ++fv_it) {

            // load color
            if(m_show_color) {

                OpenMesh::Vec3uc color = m_mesh->color(fv_it);
                glColor3ub(color[0],color[1],color[2]);

            }
            else
                glColor3ub(200,200,200);


            // load vertex
            glVertex3fv(&m_mesh->point(fv_it)[0]);

        }

    }

    glEnd();

}*/

