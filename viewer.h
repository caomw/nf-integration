#ifndef VIEWER_H
#define VIEWER_H

#include <QGLWidget>
#include <QWheelEvent>

#include "cam.h"
#include "trimesh.h"

class CViewer:public QGLWidget {

    Q_OBJECT

public:

    //! Constructor.
    explicit CViewer(const CPinholeCam<float>& cam, QWidget* parent = 0);

    //! Computes depth map for a view.
    CDenseArray<float> getDepthMap(const CRigidMotion<float,3>& viewpoint);

    //! Method stump.
    virtual void updateBoundingBox() {}

signals:

public slots:

    //! A slot that processes signals that indicate a change in view point.
    virtual void updateView(const CViewPoint<float>& viewpoint) { m_viewpoint = viewpoint; loadView(m_viewpoint.GetTransformation()); updateGL(); }

    //! Use to trigger changes of the color flag from outside.
    void onShowErrorChanged(int newstate) { m_show_color = bool(newstate);  updateGL(); }

protected:

    //! Initializes OpenGL context.
    void initializeGL();

    //! Triggers drawing of the scene.
    void paintGL();

    //! Translating wheel event into zooming.
    void wheelEvent(QWheelEvent* event);

    //! Handling of mouse clicks.
    void mousePressEvent(QMouseEvent* event);

    //! Handling of mouse drags.
    void mouseMoveEvent(QMouseEvent* event);

    //! Handling of keyboard inputs.
    void keyPressEvent(QKeyEvent* event);

    //! Method stump.
    virtual void updateClipDepth(const CRigidMotion<float,3>& viewpoint, float tolerance = 1.0) {}

protected:

    CViewPoint<float> m_viewpoint;         //!< vantage point
    CPinholeCam<float> m_cam;              //!< intrinsic camera parameters
    float m_znear;                         //!< near clipping plane
    float m_zfar;                          //!< far clipping plane
    QPoint m_last_point;                   //!< auxiliary variable to store mouse pointer locations
    CVector<float,3> m_center;             //!< center of camera rotations
    //CBoundingBox<float> m_bbox;          //!< bounding box of the scene
    bool m_show_color;                     //!< color flag

    //! Sends a view to OpenGL.
    void loadView(const CRigidMotion<float,3>& viewpoint);

    //! Updates projection matrix (e.e. if clip depths changed).
    void loadProjectionMatrix();


};


//class CTriMeshViewer:public CViewer {

//    Q_OBJECT

//public:

//    //! Constructor.
//    explicit CTriMeshViewer(const R4R::CView<float>& view, const CTriangleMesh* mesh = nullptr, QWidget* parent = nullptr);

//    //! Updates mesh and recompute its bounding box.
//    void setMesh(const CTriangleMesh* mesh);

//    /*! \brief Triggers update of bounding box without changing the pointer to the mesh.
//     *
//     * This comes in handy e.g. when new polygons are added to the mesh or the mesh is deformed.
//     */
//    void updateBoundingBox();

//    //! Assigns to each pixel in the result the handle of the face it sees.
//    CDenseArray<int> getFaceMap(const CView<float>& view);

//public slots:

//    /*! \copybrief CViewer::updateView(const R4R::CView<double>&)
//     *
//     * This needs to be overridden because we have to adjust the clip depths if the
//     * vantage point changes.
//     *
//     */
//    void updateView(const R4R::CView<float>& view);

//protected:

//    //! \copydoc CViewer::paintGL()
//    void paintGL();

//    //! Handling of mouse release event.
//    void mouseReleaseEvent(QMouseEvent* event);

//    //! Update clip depth based on the bounding box approximation of the current mesh.
//    void updateClipDepth(const CView<float>& view, float tolerance = 1.0);

//private:

//    const CTriangleMesh* m_mesh;                //!< pointer to the point cloud

//};


#endif // VIEWER_H
