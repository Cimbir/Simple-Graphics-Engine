// g++ -c .\renderer.cpp -ID:\Workspace\SFML-2.6.1\include
// g++ .\renderer.o -o Renderer -LD:\Workspace\SFML-2.6.1\lib -lgsl -lsfml-graphics -lsfml-window -lsfml-system
#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <SFML/Graphics.hpp>
#include <gsl/gsl_sf_bessel.h>
#include <conio.h>
#include <SFML/Graphics.hpp>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>

using namespace std;
using namespace sf;

#define edge pair<int,int>
#define screenPointResult pair<Vector2f,bool>

int width = 1000;
int height = 1000;

Vector3f CENTER(-10,-10,30);
double THETA = M_PI_4;
double GAMMA = -M_PI_4;

Vector3f normal;
Vector3f vertical;
Vector3f horizontal;

double moveSpeed = 0.1;
double rotateSpeed = 0.1;

double distanceFromScreen = 2;
double verticalScreenRadius = 1;
double horizontalScreenRadius = 1;

static Vector3f operator*(Vector3f a, double b){
    Vector3f res = a;
    res.x *= b;
    res.y *= b;
    res.z *= b;
    return res;
}
void updateVectors(){
    normal = Vector3f(
        cos(THETA) * cos(GAMMA),
		sin(THETA) * cos(GAMMA),
		sin(GAMMA)
    );
    vertical = Vector3f(
        cos(THETA) * cos(GAMMA - M_PI_2),
		sin(THETA) * cos(GAMMA - M_PI_2),
		sin(GAMMA - M_PI_2)
    );
    horizontal = Vector3f(
        cos(THETA - M_PI_2),
		sin(THETA - M_PI_2),
		0
    );
}
screenPointResult getPointPositionOnScreen(Vector3f p){
    Vector3f centerToPoint = p-CENTER;

    // h1 v1 n1   x   a1
    // h2 v2 n2 * y = a2 * k
    // h3 v3 n3   d   a3
    // h - horizontal   |
    // v - vertival     | A
    // n - normal       |
    // (x,y) - position on screen
    // d - distanceFromScreen
    // a - centerToPoint
    // k - positive coef

    // Invert A
    double A[] = {
        horizontal.x, vertical.x, normal.x,
        horizontal.y, vertical.y, normal.y,
        horizontal.z, vertical.z, normal.z
    };
    gsl_matrix_view AM = gsl_matrix_view_array(A, 3, 3);
    double invA[9];
    int s;
    gsl_matrix_view invAM = gsl_matrix_view_array(invA, 3, 3);
    gsl_permutation* perm = gsl_permutation_alloc(3);

    try{
        gsl_linalg_LU_decomp(&AM.matrix, perm, &s);
        gsl_linalg_LU_invert(&AM.matrix, perm, &invAM.matrix);
    }catch(exception ex){ }

    gsl_permutation_free(perm);

    // invA * a
    double a[] = { centerToPoint.x, centerToPoint.y, centerToPoint.z };
    gsl_matrix_view aM = gsl_matrix_view_array(a, 3, 1);

    double invAa[] = { 0, 0, 0 };
    gsl_matrix_view invAaM = gsl_matrix_view_array(invAa, 3, 1);

    try{
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &invAM.matrix, &aM.matrix, 0.0, &invAaM.matrix);
    }catch(exception ex){
        cout << "Here" << endl;
    }

    // finding k
    double k = distanceFromScreen / invAa[2];
    // finding (x,y) 
    Vector2f res(invAa[0] * k, invAa[1] * k);
    return {res, k >= 0};
}
Vector2f pointIfCutsScreen(Vector3f a, Vector3f b){
    // a -> b
    Vector3f u = b-a;

    double coefs[] = {
        horizontal.x, vertical.x, u.x,
        horizontal.y, vertical.y, u.y,
        horizontal.z, vertical.z, u.z
    };
    gsl_matrix_view coefsM = gsl_matrix_view_array(coefs, 3, 3);
    double res[] = {
        a.x - CENTER.x - distanceFromScreen * normal.x,
        a.y - CENTER.y - distanceFromScreen * normal.y,
        a.z - CENTER.z - distanceFromScreen * normal.z
    };
    gsl_vector_view resM = gsl_vector_view_array(res, 3);

    double x[] = {0,0,0};
    gsl_vector_view xm = gsl_vector_view_array(x, 3);

    int s;
    gsl_permutation* p = gsl_permutation_alloc(3);

    try{
        gsl_linalg_LU_decomp(&coefsM.matrix, p, &s);
        gsl_linalg_LU_solve(&coefsM.matrix, p, &resM.vector, &xm.vector);
    }catch(exception ex){
        cout << "Here" << endl;
    }

    gsl_permutation_free(p);
    return Vector2f(x[0], x[1]);
}
Vector2f screenPosToCoords(Vector2f p){
    Vector2f res;
    res.x = width/2 + p.x * (width/(2*horizontalScreenRadius));
    res.y = height/2 + p.y * (height/(2*verticalScreenRadius));
    //cout << res.x << " " << res.y << endl;
    return res;
}
void detectKeyboardInput(){
    if(Keyboard::isKeyPressed(Keyboard::A)) CENTER -= horizontal * moveSpeed;
    if(Keyboard::isKeyPressed(Keyboard::D)) CENTER += horizontal * moveSpeed;
    if(Keyboard::isKeyPressed(Keyboard::W)) CENTER -= vertical * moveSpeed;
    if(Keyboard::isKeyPressed(Keyboard::S)) CENTER += vertical * moveSpeed;
    if(Keyboard::isKeyPressed(Keyboard::R)) CENTER += normal * moveSpeed;
    if(Keyboard::isKeyPressed(Keyboard::F)) CENTER -= normal * moveSpeed;

    if(Keyboard::isKeyPressed(Keyboard::Left)) THETA += M_PI_4 * rotateSpeed;
    if(Keyboard::isKeyPressed(Keyboard::Right)) THETA -= M_PI_4 * rotateSpeed;
    if(Keyboard::isKeyPressed(Keyboard::Up)) GAMMA += M_PI_4 * rotateSpeed;
    if(Keyboard::isKeyPressed(Keyboard::Down)) GAMMA -= M_PI_4 * rotateSpeed;
}
vector<VertexArray> createObject(vector<Vector3f>& points, vector<edge>& connections){
    vector<Vector2f> pointsOnScreen(points.size());
    vector<bool> onScreen(points.size(), false);
    for(int i = 0; i < points.size(); i++) {
        screenPointResult tmp = getPointPositionOnScreen(points[i]);
        onScreen[i] = tmp.second;
        pointsOnScreen[i] = screenPosToCoords(tmp.first);
    }

    vector<VertexArray> res;
    for(int i = 0; i < connections.size(); i++){
        VertexArray line(Lines, 2);
        Vector2f screenCut = screenPosToCoords(pointIfCutsScreen(points[connections[i].first], points[connections[i].second]));
        line[0] = pointsOnScreen[connections[i].first]; 
        line[1] = pointsOnScreen[connections[i].second];
        if(!onScreen[connections[i].first]) {
            //cout << "HERE 1" << endl;
            line[0] = screenCut;
        }
        if(!onScreen[connections[i].second]) {
            //cout << "HERE 2" << endl;
            line[1] = screenCut;
        }
        res.push_back(line);
    }

    return res;
}
void display(vector<VertexArray>& lines, RenderWindow& window){
    for(int i = 0; i < lines.size(); i++){
        window.draw(lines[i]);
    }
}

int main(){
    gsl_set_error_handler_off();

            // Cube centered at origin
        // vector<Vector3f> cubeVertices = {
        //     Vector3f(5, 5, 5),
        //     Vector3f(5, -5, 5),
        //     Vector3f(5, -5, -5),
        //     Vector3f(5, 5, -5),
        //     Vector3f(-5, 5, 5),
        //     Vector3f(-5, -5, 5),
        //     Vector3f(-5, -5, -5),
        //     Vector3f(-5, 5, -5)
        // };
        // vector<edge> cubeEdges = {
        //     {0,1}, {1,2}, {2,3}, {3,0},
        //     {0,4}, {1,5}, {2,6}, {3,7},
        //     {4,5}, {5,6}, {6,7}, {7,4}
        // };
        // vector<VertexArray> cube = createObject(cubeVertices, cubeEdges);
        // display(cube, window);

        // Axis representations
        // vector<Vector3f> axis = {
        //     Vector3f(100, 0, 0),
        //     Vector3f(-100, 0, 0),
        //     Vector3f(0, 100, 0),
        //     Vector3f(0, -100, 0),
        //     Vector3f(0, 0, 100),
        //     Vector3f(0, 0, -100)
        // };
        // vector<edge> axisConnections = {
        //     {0,1}, {2,3}, {4,5}
        // };
        // vector<VertexArray> axisLines = createObject(axis, axisConnections);
        // display(axisLines, window);

        srand((int)time(0));
        int meshAmount = 50;
        vector<Vector3f> meshPoints;
        for(int i = 0-meshAmount/2; i < meshAmount-meshAmount/2; i++){
            for(int j = 0-meshAmount/2; j < meshAmount-meshAmount/2; j++){
                meshPoints.push_back(Vector3f(i*10, j*10, rand() % 20));
            }
        }
        vector<edge> meshEdges;
        for(int i = 0; i < meshAmount; i++){
            for(int j = 0; j < meshAmount; j++){
                if((i+1) < meshAmount) meshEdges.push_back({i*meshAmount+j, (i+1)*meshAmount+j});
                if((j+1) < meshAmount) meshEdges.push_back({i*meshAmount+j, i*meshAmount+(j+1)});
                if((i+1) < meshAmount && (j+1) < meshAmount) meshEdges.push_back({i*meshAmount+j, (i+1)*meshAmount+(j+1)});
            }
        }

    // create the window
    RenderWindow window(VideoMode(width, height), "My window");
    window.setFramerateLimit(30);

    // run the program as long as the window is open
    while (window.isOpen())
    {
        // check all the window's events that were triggered since the last iteration of the loop
        Event event;
        while (window.pollEvent(event))
        {
            // "close requested" event: we close the window
            if (event.type == Event::Closed)
                window.close();
        }

        // clear the window with black color
        window.clear(Color::Black);

        // draw everything here...
        // window.draw(...);

        updateVectors();

        vector<VertexArray> mesh = createObject(meshPoints, meshEdges);
        display(mesh, window);

        detectKeyboardInput();

        // end the current frame
        window.display();
    }

    return 0;
}