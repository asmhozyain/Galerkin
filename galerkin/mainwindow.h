#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QGraphicsScene>
#include <QVector>
#include <vector>
#include <QDebug>
#include <complex>
#include <cmath>
#include <math.h>
#include <QMatrix>
#include <stdio.h>
#include <QGraphicsScene>
#include <iostream>
#include <numeric>
#include "qmatrix.h"
#include "QGenericMatrix"
#include <Eigen/Dense>


QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:
    void on_Calc_clicked();

private:
    Ui::MainWindow *ui;
    void build(double w, double h, double d);
    void calculation (std::vector<std::vector<double>> second, int n);
    void draw(double w, double h, double d);
    QGraphicsScene *scene;

};
#endif // MAINWINDOW_H
