#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

#include <vector>

#include "nfield.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private:

    Ui::MainWindow *ui;

    /*
     * in file dialogue, check box to indicate whether this is a field or
     * a deflectometric image!
     *
     */
    std::vector<CNormalField<float>*> m_imgs;       //! multi-view normal field



};

#endif // MAINWINDOW_H
