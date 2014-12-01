//////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2014, Jonathan Balzer
//
// All rights reserved.
//
// This file is part of the software accompanying the following publication
//
// @INPROCEEDINGS{Balzeretal2014,
// author = {J. Balzer and D. Acevedo-Feliz and S. Soatto and S. H\"ofer and M. Hadwiger and J. Beyerer},
// title = {Cavlectometry: Towards Holistic Reconstruction of Large Mirror Objects},
// booktitle = {International Conference on 3D Vision (3DV)},
// year = {2014}}
//
// This file contains free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This source code is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this file. If not, see <http://www.gnu.org/licenses/>.
//
//////////////////////////////////////////////////////////////////////////////////

#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

#include <vector>

#include "nfield.h"
#include "viewer.h"

// forward declaration
namespace Ui {

class MainWindow;

}

class MainWindow : public QMainWindow {

    // tell QT to process this file with MOC
    Q_OBJECT

public:

    //! Constructor.
    explicit MainWindow(QWidget* parent = 0);

    //! Destructor.
    ~MainWindow();

signals:

    //! Change intrinsic camera parameters.
    void cameraChanged(const CCamera<float>& cam);

    //! Change vantage point.
    void viewpointChanged(const CViewPoint<float>& cam);

private slots:

    //! Quit program.
    void on_actionExit_triggered();

    //! Show OpenGL viewer.
    void on_action3D_Scene_triggered();

    //! Reload current mesh.
    void on_actionReload_triggered();

    //! Load normal fields.
    void on_actionOpen_View_triggered();

    //! Load initial triangular mesh.
    void on_actionLoad_Mesh_triggered();

    //! Iterate trough normal fields acquired from different vantage points.
    void on_imgSpinBox_valueChanged(int arg1);

    //! Carry out one step of the integration algorithm.
    void on_stepButton_clicked();

    //! Save deformed mesh to disk.
    void on_actionSave_Mesh_triggered();

    //! Compute and visualize local normal error distribution.
    void on_actionError_triggered();

    //! Show a visualization of the normal field.
    void showNormalFieldImage(int no);

    //! Uniform mesh refinement.
    void on_actionRefine_triggered();

    //! Save visualization of normal field to disk.
    void on_actionSave_Image_triggered();

    //! Visualize the region on the surface visible by the current view.
    void on_actionVisibility_triggered();

    //! Close the current view.
    void on_actionClose_Current_triggered();

    //! Apply one step of Laplacian smoothinh to mesh.
    void on_actionSmooth_triggered();

    //! Compute the views in which a face is visible.
    void computeVisibility(OpenMesh::FPropHandleT<std::set<uint> > views);

    //! Compute normal residual.
    double computeResidual(CDenseVector<double>& r);

    //! Delete parts of the mesh for which the mask is equal to zero.
    void on_actionCrop_triggered();

private:

    Ui::MainWindow* ui;                                //! the GUI
    std::vector<CNormalField<float> > m_imgs;          //! multi-view normal field
    CTriangleMesh m_mesh;                              //!< mesh to be optimizedattached to faces of mesh
    QString m_mesh_filename;                           //!< location where mesh is stored on disk
    CTriMeshViewer* m_viewer;                          //! OpenGL widget

};

#endif // MAINWINDOW_H
