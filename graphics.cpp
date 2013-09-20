#include "graphics.h"

#include <cmath>
#include <iomanip>
#include <QImage>
#include <QString>
#include <QImageWriter>
#include <QDateTime>

Graphics::Graphics(Diffusion2DpAR *_mesh, tFloat *_X, QString _str, QWidget *parent) :
    QGLWidget(parent)
{
    this->showMaximized();

    mesh = _mesh;
    X = _X;
    this->setWindowTitle(_str);

    xmin = -1.2;
    ymin = -1.2;
    xmax = +1.2;
    ymax = +1.2;

}


void Graphics::initializeGL()
{
    glShadeModel(GL_SMOOTH);

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glDepthFunc(GL_LEQUAL);

    glEnable(GL_POINT_SMOOTH);
    glHint(GL_POINT_SMOOTH_HINT,GL_NICEST);

    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);

    glClearColor(1.0f, 1.0f, 1.0f, 0.0f);
    glClearDepth(1.0f);

    glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);

}

void Graphics::resizeGL(int width, int height)
{
    double width_ = static_cast<double>(this->width());
    double height_ = static_cast<double>(this->height());

    if (width_ > height_)
    {
        height_ = height_?height_:1;
        double correction = 0.5 * ( width_/ height_ * (ymax-ymin) - (xmax-xmin) );
        xmin   -= correction;
        xmax +=correction;
    }
    else
    {
        width_ = width_?width_:1;
        double correction = 0.5 * ( height_ / width_ * (xmax-xmin) - (ymax-ymin) );
        ymin   -= correction;
        ymax  += correction;
    }

    glViewport( 0, 0, (GLint)width, (GLint)height );

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(xmin, xmax, ymin, ymax, 1.0, -1.0);
    //gluPerspective(60, (float)width/height, 0.1, 50000);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

}


void Graphics::wheelEvent(QWheelEvent *event)
{
    double zoom_factor = -event->delta()/120*0.2;

    double X_opengl = event->x()/static_cast<double>(this->width())*(xmax - xmin)+xmin;
    double Y_opengl  = (this->height()-event->y())/static_cast<double>(this->height())*(ymax - ymin)+ymin;

    xmin -= (X_opengl-xmin)*zoom_factor;
    xmax += (xmax-X_opengl)*zoom_factor;

    ymin -= (Y_opengl-ymin)*zoom_factor;
    ymax += (ymax-Y_opengl)*zoom_factor;


    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    glOrtho(xmin, xmax, ymin, ymax, 1.0, -1.0);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    this->repaint();

}

void Graphics::mousePressEvent(QMouseEvent *event)
{
    isMousePress = true;
    panX = event->x();
    panY = event->y();
}


void Graphics::mouseReleaseEvent(QMouseEvent *event)
{
    isMousePress = false;
}

void Graphics::mouseDoubleClickEvent(QMouseEvent *event)
{
    QDateTime now = QDateTime::currentDateTime();

    QString filename = QString("imc-diffusion2d-snapshot-")
            + now.toString("yyyyMMddhhmmsszzz") + QString(".png");
    this->updateGL();
    this->grabFrameBuffer(true).save(filename, "PNG", 100);
}


void Graphics::mouseMoveEvent(QMouseEvent *event)
{
    if(event->buttons() == Qt::LeftButton)
    {
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();

        double X_opengl = (-event->x()+panX)/static_cast<double>(this->width())*(xmax - xmin);
        double Y_opengl  = (event->y()-panY)/static_cast<double>(this->height())*(ymax - ymin);

        xmax += X_opengl;
        xmin += X_opengl;

        ymax += Y_opengl;
        ymin += Y_opengl;

        panX = event->x();
        panY = event->y();

        glOrtho(xmin, xmax, ymin, ymax, 1.0, -1.0);

        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();

        this->repaint();

    }

    updateGL();

}

void Graphics::paintGL()

{
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();


    glColor4d(0.0, 0.0, 0.0, 1.0);
    glBegin(GL_QUADS);
    {
        glVertex3d(0.0, 0.0, 0.0);
        glVertex3d(QtoD(mesh->lx), 0.0, 0.0);
        glVertex3d(QtoD(mesh->lx), QtoD(mesh->ly), 0.0);
        glVertex3d(0.0, QtoD(mesh->ly), 0.0);
    }
    glEnd();


    // Desenhar os Nós
    glPointSize(5.0f);
    glColor4d(1.0, 1.0, 1.0, 1.0);
    glBegin(GL_POINTS);
    for(int i=0; i<mesh->nnodes; i++){
        if(mesh->nodes[i].boundary == NULL)
            glColor4d(1.0, 1.0, 1.0, 1.0);
        else
            glColor4d(0.0, 0.0, 0.0, 1.0);
        glVertex2d(QtoD(mesh->nodes[i].x), QtoD(mesh->nodes[i].y));
    }
    glEnd();

    //    // Desenhar os Elementos
    //    glColor4d(1.0, 1.0, 1.0, 1.0);

    //    glBegin(GL_LINES);
    //    for(int i=0; i<mesh->nx; i++){
    //        glVertex2d(QtoD(i*mesh->hx), QtoD(0.0));
    //        glVertex2d(QtoD(i*mesh->hx), QtoD(mesh->ly));
    //    }
    //    for(int i=0; i<mesh->ny; i++){
    //        glVertex2d(QtoD(0.0), QtoD(i*mesh->hy));
    //        glVertex2d(QtoD(mesh->lx), QtoD(i*mesh->hy));
    //    }
    //    glEnd();


    //    // Desenhar os Nós
    //    glPointSize(10.0f);
    //    glColor4d(0.0, 1.0, 1.0, 1.0);
    //    glBegin(GL_POINTS);

    //    for(int i=0; i<mesh->nnodes; i++){
    //        if(mesh->nodes[i].nodes[0] == NULL ||
    //                mesh->nodes[i].nodes[1] == NULL ||
    //                mesh->nodes[i].nodes[2] == NULL ||
    //                mesh->nodes[i].nodes[3] == NULL
    //                )

    //            glVertex2d(QtoD(mesh->nodes[i].x), QtoD(mesh->nodes[i].y));
    //    }
    //    glEnd();



    // Desenha os resultados
    double T0, T1, T2, T3, T4, Tn, R, G, B;
    getMaxMin(X, static_cast<int>(mesh->nnodes), T4, T0);

    T2 = (T0+T4)/2.;
    T1 = (T0+T2)/2.;
    T3 = (T2+T4)/2.;

    int k[4];

    //    for(int i=0; i<mesh->nx*(mesh->ny-1); i++){

    //        k[0] = i;
    //        k[1] = i+1;
    //        k[2] = i+mesh->nx+1;
    //        k[3] = i+mesh->nx;

    //        glBegin(GL_QUADS);
    //        for(int p = 0; p<4; p++){
    //            Tn = X[k[p]];
    //            R = Tn<T2?  0. : (Tn>T3? 1. : (Tn-T2)/(T3-T2));
    //            B = Tn>T2?  0. : (Tn<T1? 1. : (T2-Tn)/(T2-T1));
    //            G = Tn<T1? (Tn-T0)/(T1-T0) : Tn>T3 ? (T4-Tn)/(T4-T3) : 1.;
    //            glColor4d(R,G,B,0.8);

    //            glVertex2d(mesh->nodes[k[p]].x, mesh->nodes[k[p]].y);
    //        }
    //        glEnd();

    //        if(i+1 && !((i+2)%mesh->nx))
    //            i+=1;
    //    }


    //    for(int i=0; i<mesh->nelements; i++){

    //        if(mesh->elements[i].index == -1) continue;

    //        k[0] = mesh->elements[i].nodes[0]->index;
    //        k[1] = mesh->elements[i].nodes[1]->index;
    //        k[2] = mesh->elements[i].nodes[2]->index;
    //        k[3] = mesh->elements[i].nodes[3]->index;

    //        glBegin(GL_QUADS);
    //        for(int p = 0; p<4; p++){
    //            Tn = X[k[p]];
    //            R = Tn<T2?  0. : (Tn>T3? 1. : (Tn-T2)/(T3-T2));
    //            B = Tn>T2?  0. : (Tn<T1? 1. : (T2-Tn)/(T2-T1));
    //            G = Tn<T1? (Tn-T0)/(T1-T0) : Tn>T3 ? (T4-Tn)/(T4-T3) : 1.;
    //            glColor4d(R,G,B,0.8);

    //            glVertex2d(mesh->nodes[k[p]].x, mesh->nodes[k[p]].y);
    //        }
    //        glEnd();

    //    }

    int pk[100];
    int count;
    Node2D *I, *F, *V;

    for(int i=0; i<mesh->nelements; i++){


        if(mesh->elements[i].nodes[0] == NULL||
                mesh->elements[i].nodes[1] == NULL ||
                mesh->elements[i].nodes[2] == NULL ||
                mesh->elements[i].nodes[3] == NULL)
            continue;

        if(mesh->elements[i].index == -1) continue;

        count = 0;
        I = mesh->elements[i].nodes[0];
        F = mesh->elements[i].nodes[1];
        pk[count++] = I->index;
        do{
            I = I->nodes[1];
            pk[count++] = I->index;
        }while(I!=F);

        I = mesh->elements[i].nodes[1];
        F = mesh->elements[i].nodes[2];
        //pk[count++] = I->index;
        do{
            I = I->nodes[2];
            pk[count++] = I->index;
        }while(I!=F);

        I = mesh->elements[i].nodes[2];
        F = mesh->elements[i].nodes[3];
        //pk[count++] = I->index;
        do{
            I = I->nodes[3];
            pk[count++] = I->index;
        }while(I!=F);

        I = mesh->elements[i].nodes[3];
        F = mesh->elements[i].nodes[0];
        //pk[count++] = I->index;
        do{
            I = I->nodes[0];
            pk[count++] = I->index;
        }while(I!=F);


        //        for(int j=1; j<count-1; j++){
        //            k[0] = pk[0];
        //            k[1] = pk[j];
        //            k[2] = pk[j+1];
        //            glBegin(GL_TRIANGLES);
        //            for(int p = 0; p<3; p++){
        //                Tn = X[k[p]];
        //                R = Tn<T2?  0. : (Tn>T3? 1. : (Tn-T2)/(T3-T2));
        //                B = Tn>T2?  0. : (Tn<T1? 1. : (T2-Tn)/(T2-T1));
        //                G = Tn<T1? (Tn-T0)/(T1-T0) : Tn>T3 ? (T4-Tn)/(T4-T3) : 1.;
        //                glColor4d(R,G,B,0.8);

        //                glVertex2d(mesh->nodes[k[p]].x, mesh->nodes[k[p]].y);
        //            }
        //            glEnd();
        //        }

        glBegin(GL_TRIANGLE_FAN);
        for(int j=0; j<count; j++){
            k[0] = pk[0];
            k[1] = pk[j];
            k[2] = pk[j+1];

            //for(int p = 0; p<3; p++){
            Tn = X[pk[j]];
            R = Tn<T2?  0. : (Tn>T3? 1. : (Tn-T2)/(T3-T2));
            B = Tn>T2?  0. : (Tn<T1? 1. : (T2-Tn)/(T2-T1));
            G = Tn<T1? (Tn-T0)/(T1-T0) : Tn>T3 ? (T4-Tn)/(T4-T3) : 1.;
            glColor4d(R,G,B,0.8);

            glVertex2d(mesh->nodes[pk[j]].x, mesh->nodes[pk[j]].y);
            //}

        }
        glEnd();

    }



    //    for(tInteger i=0; i<30; i++)


    //        if(mesh->nodes[i].nodes[1] != NULL
    //                && mesh->nodes[i].nodes[1]->nodes[2] !=NULL  && mesh->nodes[i].nodes[2] != NULL){

    //        k[0] = mesh->nodes[i].index;
    //        k[1] = mesh->nodes[i].nodes[1]->index;
    //        k[2] = mesh->nodes[i].nodes[1]->nodes[2]->index;
    //        k[3] = mesh->nodes[i].nodes[2]->index;

    //        glBegin(GL_QUADS);
    //        for(int p = 0; p<4; p++){
    //            Tn = X[k[p]];
    //            R = Tn<T2?  0. : (Tn>T3? 1. : (Tn-T2)/(T3-T2));
    //            B = Tn>T2?  0. : (Tn<T1? 1. : (T2-Tn)/(T2-T1));
    //            G = Tn<T1? (Tn-T0)/(T1-T0) : Tn>T3 ? (T4-Tn)/(T4-T3) : 1.;
    //            glColor4d(R,G,B,0.8);

    //            glVertex2d(mesh->nodes[k[p]].x, mesh->nodes[k[p]].y);
    //        }
    //        glEnd();

    //    }

    glColor4d(1.0, 1.0, 1.0, 0.5);
    for(tInteger i=0; i<mesh->nnodes; i++)
    {
        glBegin(GL_LINES);
        for(tInteger j=0; j<4;j++)
            if(mesh->nodes[i].nodes[j] != NULL)
            {
                glVertex2d(mesh->nodes[i].x, mesh->nodes[i].y);
                glVertex2d(mesh->nodes[i].nodes[j]->x, mesh->nodes[i].nodes[j]->y);
            }

        glEnd();
    }




    this->drawLegend(T0, T1, T2, T3, T4);

    for(tInteger i=0; i<mesh->nnodes; i++)
    {
        QString st0 = QString("%1").arg(static_cast<int>(mesh->nodes[i].index));
        glColor3f(0.,0.,0.);
        QFont font10("Ubuntu", 8, QFont::Cursive);
        this->renderText(QtoD(mesh->nodes[i].x), QtoD(mesh->nodes[i].y) ,0.f,st0,font10);


    }

}


void Graphics::getMaxMin(tFloat *vector, int size, double &max, double &min)
{
    max = vector[0];
    min = vector[0];

    for(int i=1; i<size; i++){
        max = vector[i]>max ? vector[i] : max;
        min = vector[i]<min ? vector[i] : min;
    }
}


void Graphics::drawLegend(double T0, double T1, double T2, double T3, double T4)
{
    //glPushMatrix();
    //glTranslated(xmin + 0.9*(xmax-xmin), ymin + 0.15*(ymax-ymin), 0.0);
    glTranslated(1.2, 0.0, 0.0);

    //glPushMatrix();
    //glScaled(0.2*(xmax-xmin), 0.7*(ymax-ymin), 1.0);

    int pn = 10;

    QString st0 = QString("%1").arg(T0, 0, 'E', pn);
    QString st1 = QString("%1").arg(T1, 0, 'E', pn);
    QString st2 = QString("%1").arg(T2, 0, 'E', pn);
    QString st3 = QString("%1").arg(T3, 0, 'E', pn);
    QString st4 = QString("%1").arg(T4, 0, 'E', pn);

    glColor3f(0.,0.,0.);
    QFont font10("Ubuntu", 12, QFont::Cursive);



    this->renderText(0.f,0.00f,0.f,st0,font10);
    this->renderText(0.f,0.245f,0.f,st1,font10);
    this->renderText(0.f,0.49f,0.f,st2,font10);
    this->renderText(0.f,0.74f,0.f,st3,font10);
    this->renderText(0.f,0.98f,0.f,st4,font10);

    //glLoadIdentity();
    glTranslated(-0.05, 0.0, 0.0);

    T0 = 0.f;
    T4 = 1.f;
    T2 = (T0+T4)/2.;
    T1 = (T0+T2)/2.;
    T3 = (T2+T4)/2.;

    glColor3f(0.,1.,0.);

    double Tn,R,B,G;

    double dx = 0.01;
    for(int i=0; i<100; i++)
    {
        glBegin(GL_QUADS);

        Tn = i*dx;
        R = Tn<T2?  0. : (Tn>T3? 1. : (Tn-T2)/(T3-T2));
        B = Tn>T2?  0. : (Tn<T1? 1. : (T2-Tn)/(T2-T1));
        G = Tn<T1? (Tn-T0)/(T1-T0) : Tn>T3 ? (T4-Tn)/(T4-T3) : 1.;
        glColor3f(R,G,B);
        glVertex3f(0.f,i*dx,0.f);
        glVertex3f(0.03f,i*dx,0.f);

        Tn = (i+1)*dx;
        R = Tn<T2?  0. : (Tn>T3? 1. : (Tn-T2)/(T3-T2));
        B = Tn>T2?  0. : (Tn<T1? 1. : (T2-Tn)/(T2-T1));
        G = Tn<T1? (Tn-T0)/(T1-T0) : Tn>T3 ? (T4-Tn)/(T4-T3) : 1.;
        glColor3f(R,G,B);
        glVertex3f(0.03f,(i+1)*dx,0.f);
        glVertex3f(0.f,(i+1)*dx,0.f);

        glEnd();
    }
    //glPopMatrix();
}

