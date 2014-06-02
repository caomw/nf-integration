#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

#include <vector>

#include "nfield.h"
#include "viewer.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

signals:
    void cameraChanged(const CCamera<float>& cam);
    void viewpointChanged(const CViewPoint<float>& cam);

private slots:
    void on_actionExit_triggered();

    void on_action3D_Scene_triggered();

    void on_actionReload_triggered();

    void on_actionOpen_View_triggered();

    void on_actionLoad_Mesh_triggered();

private:

    Ui::MainWindow* ui;                                //! the GUI
    std::vector<CNormalField<float> > m_imgs;          //! multi-view normal field
    CTriangleMesh m_mesh;                              //!< mesh to be optimized
    QString m_mesh_filename;                           //!< location where mesh is stored on disk
    OpenMesh::FPropHandleT<vec3f> m_target_normal;     //!< holds the target normal field for integration
    CTriMeshViewer* m_viewer;                          //! OpenGL widget

};

#endif // MAINWINDOW_H
