#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <thread>
#include <QDateTime>

constexpr double epsilon() {return 8.85418782 * pow(10, -12);}

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    scene = new QGraphicsScene();
    ui->graphicsView->setScene(scene);
    ui->graphicsView->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    ui->graphicsView->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
}

MainWindow::~MainWindow()
{
    delete ui;
}


void MainWindow::on_Calc_clicked()
{
    ui->Final->clear();

    QDateTime start = QDateTime::currentDateTime();

    double w = ui->Wparam->toPlainText().toDouble();
    double h = ui->Hparam->toPlainText().toDouble();
    double d = ui->Dparam->toPlainText().toDouble();

    std::thread buildT([&](){
        build(w, h, d);
    });
    //build(w, h, d); // построение структуры (координаты)




    draw(w, h, d); // отрисовка структуры в сцене

    buildT.join();

    QDateTime finish = QDateTime::currentDateTime();

    int secs = finish.secsTo(start);
    start.addSecs(secs);
    int msecs = finish.time().msecsTo(start.time());

    int msecs_duration = secs * 1000 + msecs;
    qDebug() << msecs_duration;



}

void MainWindow::build(double w, double h, double d)
{

    double x, y, len;

    int Segment = ui->Nseg->toPlainText().toInt();
    std::vector<double> first1;
    std::vector<double> first2;
    std::vector<std::vector<double>> second1;
    std::vector<std::vector<double>> second2;
    int n = 4;



    std::thread t1 ([&](){
        for (int i = 0; i < 2 ; i++) {
            for (int j = 0; j < Segment; j++) {
                first1.clear();


                    if (i == 0){
                        len = w / Segment;
                        x = (len * j) + d + len / 2;
                        y = h;

                    }

                    if (i == 1){
                        len = d / Segment;
                        x = (len * j) + len / 2;
                        y = h;
                    }

                    if (i == 2){
                        len = d / Segment;
                        x = (len * j) + d + w + len / 2;
                        y = h;
                    }

                    if (i == 3){
                        len = (w + d * 2) / Segment;
                        x = (len * j) + len / 2;
                        y = 0;
                    }

                    first1.push_back(x);
                    first1.push_back(y);
                    first1.push_back(len);
                    second1.push_back(first1);
            }

        }
    });

    std::thread t2 ([&](){
        for (int i = 2; i < n ; i++) {
            for (int j = 0; j < Segment; j++) {
                first2.clear();


                    if (i == 0){
                        len = w / Segment;
                        x = (len * j) + d + len / 2;
                        y = h;

                    }

                    if (i == 1){
                        len = d / Segment;
                        x = (len * j) + len / 2;
                        y = h;
                    }

                    if (i == 2){
                        len = d / Segment;
                        x = (len * j) + d + w + len / 2;
                        y = h;
                    }

                    if (i == 3){
                        len = (w + d * 2) / Segment;
                        x = (len * j) + len / 2;
                        y = 0;
                    }

                    first2.push_back(x);
                    first2.push_back(y);
                    first2.push_back(len);
                    second2.push_back(first2);
            }

        }
    });

    t1.join();
    t2.join();

        second1.insert(second1.end(), second2.begin(), second2.end());

    //qDebug() << second;
    calculation(second1, n * Segment);
}

void MainWindow::calculation(std::vector<std::vector<double> > second, int n)
{
    double er = ui->Erparam->toPlainText().toDouble();
    int Segment = ui->Nseg->toPlainText().toInt();

    using namespace std;

    Eigen::MatrixXcd Zz(n,n);
    Eigen::MatrixXcd Vv(n,1);
    Eigen::MatrixXd Znew (n, n);
    Eigen::MatrixXcd Aa(n,1);

    double x1, x2, y1, y2, a, x, y,len;

    for (int i = 0; i < n; i++) {

        a = second[i][2] / 2; // Область интегрирования элемента (относительно которого считаем)

        if (i < Segment){

            Vv(i,0) = 1 * second[i][2];
        }
        if (i >= Segment){

            Vv(i,0) = 0;
        }


        for (int j = 0; j < n; j++){

            x = second[j][0]; // х координата сегмента
            y = second[j][1]; // y координата сегмента
            len = second[j][2]; // длина сегмента

            x1 = x - (len / 2) - second[i][0];
            x2 = x1 + len;
            y1 = second[i][1] - y + pow(10, -12);
            y2 = y1;

            complex <double> z1(x1, y1);
            complex <double> z2(x2, y2);
            complex <double> uz = (z2 - z1) / (abs(z2 - z1));

            complex <double> Kv1 = (pow((z2 - a), 2) / 2.0) * (log(z2 - a) - 1.5);
            complex <double> Kv2 = (pow((z1 - a), 2) / 2.0) * (log(z1 - a) - 1.5);
            complex <double> Kv3 = (pow((z2 - -a), 2) / 2.0) * (log(z2 - -a) - 1.5);
            complex <double> Kv4 = (pow((z1 - -a), 2) / 2.0) * (log(z1 - -a) - 1.5);

            complex <double> solve = uz * ((Kv1 - Kv2) - (Kv3 - Kv4));


            Zz(i, j) = ((1 / (2 * M_PI * epsilon())) * solve.real()) / second[i][2];



        }


    }

    //Вычитание земли из элемента
    for(int i = 0; i < n - 1; i++){

        for (int j = 0; j < n ; j++) {

            Zz(i, j) = Zz(i, j) - Zz(n - 1,j);

        }
    }

    double temp = 0;

    // Последний элемент равен сумме элемента умноженной на длину сегмента * -1

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                temp = Zz(i, j).real() * second[j][2] + temp;
            }
            Zz(i, n - 1) = -temp;
            temp = 0;

        }





    // EIGEN SLAE

    Aa = Zz.partialPivLu().solve(Vv);
    std::vector<std::complex<double>> Av (Aa.data(), Aa.data()+Aa.rows() * Aa.cols());
    qDebug() << "Eigen solve";
    double SumE = 0;
    for (int i = 0; i < Segment; i++) {
        SumE = SumE + Av[i].real() * er;
    }

    //for(int i=0; i < Av.size(); i++) qDebug() << Av[i].real();


    ui->Final->insertPlainText(QString::number(SumE));



    //Решение СЛАУ

    //прямой ход
//    double v;
//    for (int k = 0, i, j, im; k < n - 1; k++) {
//        im = k;
//        for (i = k + 1; i < n; i++) {
//            if(fabs(Zz(im, k)) < fabs(Zz(i, k))){
//                im = i;
//            }
//        }
//        if (im != k){
//            for (j = 0; j < n; j++) {
//                v = Zz (im, j).real();
//                Zz(im, j) = Zz(k, j);
//                Zz(k, j) = v;
//            }
//            v = Vv(im).real();
//            Vv(im) = Vv(k);
//            Vv(k) = v;
//        }
//        for (i = k + 1; i < n; i++) {
//            v = 1.0 * Zz(i, k).real() / Zz(k, k).real();
//            Zz(i, k) = 0;
//            Vv(i) = Vv(i) - v * Vv(k);
//            if(v != 0){
//                for (j = k + 1; j < n; j++) {
//                    Zz(i, j) = Zz(i, j) - v * Zz(k, j);
//                }
//            }
//        }
//    }

//    //обратный ход
//    double s = 0;
//    std::vector<double> A(n);
//    Aa(n - 1) = 1.0 * Vv(n - 1) / Zz(n - 1, n - 1);
//    for (int i = n - 2, j; 0 <= i; i--) {
//        s = 0;
//        for (j = i + 1; j < n; j++) {
//            s = s + Zz(i, j).real() * Aa(j).real();
//        }
//        Aa(i) = 1.0 * (Vv(i) - s) / Zz(i, i);
//    }

//    double sum = 0;

//    for (int i = 0; i < Segment; i++) {
//        sum = sum + Aa(i).real() * er;
//    }
//    qDebug() << "GaussSolve";
//    qDebug() << sum;


//ui->Final->insertPlainText(QString::number(sum));

}

void MainWindow::draw(double w, double h, double d)
{
    double zoom;
    scene->clear();
    scene->setSceneRect(0, 0, 900, 500);
    scene->clearFocus ();

    QPen redpen(Qt::red);
    redpen.setWidth(5);
    redpen.setCapStyle(Qt::RoundCap);

    QPen rectpen(Qt::black);
    rectpen.setWidth(2);
    rectpen.setCapStyle(Qt::RoundCap);

    QBrush blackBrush (Qt::gray);

    zoom = ((900 / 100) * 90) / (w + d * 2); // рисунок занимает 70% от сцены
    if (h * zoom >= (500 / 100) * 90){
        zoom = ((500 / 100) * 90) / h;
    }

    scene->addRect((900 / 2) - zoom * (w + d * 2) / 2, 500 / 2 - zoom * 0.5 * h, zoom * (w + d * 2), zoom * h, rectpen, blackBrush); // диэлектрическая подложка
    scene->addLine((900 / 2) - zoom * w / 2, 500 / 2 - zoom * 0.5 * h, (900 / 2) + zoom * w / 2, 500 / 2 - zoom * 0.5 * h, redpen); // проводник
}




