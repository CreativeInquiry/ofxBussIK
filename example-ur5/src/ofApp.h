#pragma once

#include "ofMain.h"
#include "ofxBussIK.h"

#define MAX_NUM_NODE	1000
#define MAX_NUM_THETA	1000
#define MAX_NUM_EFFECT	100
#define RADIAN(X)	((X)*RadiansToDegrees)

enum Method {
    IK_JACOB_TRANS=0,
    IK_PURE_PSEUDO,
    IK_DLS,
    IK_SDLS ,
    IK_DLS_SVD
};

class ofApp : public ofBaseApp{
public:

    
    void setup();
    void update();
    void draw();
    
    void exit();
    void keyPressed(int key);
    void keyReleased(int key);
    void mouseMoved(int x, int y );
    void mouseDragged(int x, int y, int button);
    void mousePressed(int x, int y, int button);
    void mouseReleased(int x, int y, int button);
    void mouseEntered(int x, int y);
    void mouseExited(int x, int y);
    void windowResized(int w, int h);
    void dragEvent(ofDragInfo dragInfo);
    void gotMessage(ofMessage msg);
    
    
    void BuildKukaIIWAShape();
    void drawTree(Node* node, ofNode& tr);
    void getNode(const Node* node, ofNode& act);
    void Reset(Tree &tree, Jacobian* m_ikJacobian);
    void UpdateTargets( double T2, Tree & treeY);
    void DoUpdateStep(double Tstep, Tree & treeY, Jacobian *jacob, int ikMethod);
    
    double T = 0;
    VectorR3 targetaa[MAX_NUM_EFFECT];

    ofEasyCam cam;

    ofNode target;
    int m_ikMethod;
    Tree m_ikTree;
    vector<Node*> m_ikNodes;
    Jacobian* m_ikJacobian;
    
    vector<int> m_movingInstances;
    int m_targetInstance;
    enum
    {
        numCubesX = 20,
        numCubesY = 20
    };
};
