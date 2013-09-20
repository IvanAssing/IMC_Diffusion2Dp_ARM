#ifndef GRAPHICS_H
#define GRAPHICS_H

//#include "ui_graphics.h"

#include <QGLWidget>
#include <QMouseEvent>
#include <QCursor>

#include "diffusion2dp.h"
#include "diffusion2dpar.h"

class Graphics : public QGLWidget
{
        Q_OBJECT
    public:
        explicit Graphics(Diffusion2DpAR *mesh, tFloat *X, QString str, QWidget *parent = 0);

        double xmax, xmin, ymax, ymin, panX, panY;
        bool isMousePress;

        Diffusion2DpAR *mesh;
        tFloat *X;

        void drawLegend(double T0, double T1, double T2, double T3, double T4);
        void getMaxMin(tFloat *vector, int size, double &max, double &min);

    signals:

    public slots:

        void initializeGL();
        void resizeGL(int width, int height);
        void paintGL();

        void wheelEvent(QWheelEvent *event);
        void mousePressEvent(QMouseEvent *event);
        void mouseMoveEvent(QMouseEvent *event);
        void mouseReleaseEvent(QMouseEvent *event);
        void mouseDoubleClickEvent(QMouseEvent *event);

};


#endif // GRAPHICS_H
