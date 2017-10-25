#include "ofApp.h"

// Make slowdown factor larger to make the simulation take larger, less frequent steps
// Make the constant factor in Tstep larger to make time pass more quickly
//const int SlowdownFactor = 40;
const int SlowdownFactor = 0;		// Make higher to take larger steps less frequently
const int SleepsPerStep=SlowdownFactor;
int SleepCounter=0;
//const double Tstep = 0.0005*(double)SlowdownFactor;		// Time step


int	AxesList;		/* list to hold the axes		*/
int	AxesOn;			/* ON or OFF				*/

float	Scale, Scale2;		/* scaling factors			*/



int JointLimitsOn;
int RestPositionOn;
int UseJacobianTargets1;


int numIteration = 1;
double error = 0.0;
double errorDLS = 0.0;
double errorSDLS = 0.0;
double sumError = 0.0;
double sumErrorDLS = 0.0;
double sumErrorSDLS = 0.0;

#ifdef _DYNAMIC
bool initMaxDist = true;
extern double Excess[];
extern double dsnorm[];
#endif

//--------------------------------------------------------------
void ofApp::setup(){
    BuildKukaIIWAShape();
    m_ikJacobian = new Jacobian(&m_ikTree);
    Reset(m_ikTree,m_ikJacobian);
    cam.setNearClip(0.001);
    cam.lookAt(ofVec3f(0, 0, 0), ofVec3f(0, 0, 1));
    //    cam.setFarClip(10);
    
    m_ikMethod = IK_JACOB_TRANS;
    ofBackground(0, 0, 0);
}

//--------------------------------------------------------------
void ofApp::update(){
   }

//--------------------------------------------------------------
void ofApp::draw(){
    DoUpdateStep(ofGetLastFrameTime(), m_ikTree, m_ikJacobian, m_ikMethod);

    cam.begin();
    ofNode act;
    getNode(m_ikTree.GetRoot(), act);
    ofSetColor(255, 255, 255);
    drawTree(m_ikTree.GetRoot(), act);
    ofSetColor(255, 0, 255);
    target.draw();
    
    cam.end();
}

void ofApp::Reset(Tree &tree, Jacobian* m_ikJacobian)
{
    
    AxesOn = false;
    
    Scale  = 1.0;
    Scale2 = 0.0;		/* because add 1. to it in Display()	*/
    
    JointLimitsOn = true;
    RestPositionOn = false;
    UseJacobianTargets1 = false;
    
    
    tree.Init();
    tree.Compute();
    m_ikJacobian->Reset();
    
}

// Update target positions

void ofApp::UpdateTargets( double T2, Tree & treeY) {
    double T3 = T2 / 5.;
    targetaa[0].Set(0.6*cos(T3), 0.9*sin(T3/2), 0.15);
    target.setPosition(0.6*cos(T3),  0.9*sin(T3/2), 0.15);
}


// Does a single update (on one kind of tree)
void ofApp::DoUpdateStep(double Tstep, Tree & treeY, Jacobian *jacob, int ikMethod) {
    if ( SleepCounter==0 ) {
        T += Tstep;
        UpdateTargets( T , treeY);
    }
    
    if ( UseJacobianTargets1 ) {
        jacob->SetJtargetActive();
    }
    else {
        jacob->SetJendActive();
    }
    jacob->ComputeJacobian(targetaa);						// Set up Jacobian and deltaS vectors
    
    // Calculate the change in theta values
    switch (ikMethod) {
        case IK_JACOB_TRANS:
            jacob->CalcDeltaThetasTranspose();		// Jacobian transpose method
            break;
        case IK_DLS:
            jacob->CalcDeltaThetasDLS();			// Damped least squares method
            break;
        case IK_DLS_SVD:
            jacob->CalcDeltaThetasDLSwithSVD();
            break;
        case IK_PURE_PSEUDO:
            jacob->CalcDeltaThetasPseudoinverse();	// Pure pseudoinverse method
            break;
        case IK_SDLS:
            jacob->CalcDeltaThetasSDLS();			// Selectively damped least squares method
            break;
        default:
            jacob->ZeroDeltaThetas();
            break;
    }
    
    if ( SleepCounter==0 ) {
        jacob->UpdateThetas();							// Apply the change in the theta values
        jacob->UpdatedSClampValue(targetaa);
        SleepCounter = SleepsPerStep;
    }
    else { 
        SleepCounter--;
    }
}


void ofApp::exit()
{
    delete m_ikJacobian;
    m_ikJacobian = 0;
}


void ofApp::getNode(const Node* node, ofNode& act)
{
    ofVec3f axis(node->v.x, node->v.y, node->v.z);
    ofQuaternion rot(0, 0, 0, 1);
    if (axis.length())
    {
        rot = ofQuaternion (ofRadToDeg(node->GetTheta()), axis);
    }
    act.setOrientation(rot);
    act.setPosition(node->s.x, node->s.y, node->s.z);
    ofDrawLine(act.getPosition(), act.getPosition() + axis);
}
void ofApp::drawTree(Node* node, ofNode& tr)
{
    if (node != 0) {
        getNode(node, tr);
 
        tr.draw();
        if (node->left) {
            ofNode act;
            drawTree(node->left, act);		// Draw tree of children recursively
        }
        ofFill();
    }
}

void ofApp::BuildKukaIIWAShape()
{
    const VectorR3& unitx = VectorR3::UnitX;
    const VectorR3& unity = VectorR3::UnitY;
    const VectorR3& unitz = VectorR3::UnitZ;
    const VectorR3 unit1(sqrt(14.0) / 8.0, 1.0 / 8.0, 7.0 / 8.0);
    const VectorR3& zero = VectorR3::Zero;
    
    float minTheta = -PI;
    float maxTheta = PI;
    
    m_ikNodes.resize(7);//6DOF+additional endeffector
    
    
    m_ikNodes[0] = new Node(VectorR3(0, 0, 0), unitz, 0.08, JOINT, -2*PI, 2*PI, ofDegToRad(0.));
    m_ikTree.InsertRoot(m_ikNodes[0]);
    
    m_ikNodes[1] = new Node(VectorR3(0, -0.072238, 0.083204), -unity, 0.08, JOINT, minTheta, maxTheta, ofDegToRad(-90.));
    m_ikTree.InsertLeftChild(m_ikNodes[0], m_ikNodes[1]);
    
    m_ikNodes[2] = new Node(VectorR3(0,-0.077537,0.51141), -unity, 0.08, JOINT, minTheta, maxTheta, ofDegToRad(0.));
    m_ikTree.InsertLeftChild(m_ikNodes[1], m_ikNodes[2]);
    
    m_ikNodes[3] = new Node(VectorR3(0, -0.070608, 0.903192), -unity, 0.08, JOINT, minTheta, maxTheta, ofDegToRad(-90.));
    m_ikTree.InsertLeftChild(m_ikNodes[2], m_ikNodes[3]);
    
    m_ikNodes[4] = new Node(VectorR3(0, -0.117242, 0.950973), unitz, 0.08, JOINT, minTheta, maxTheta, ofDegToRad(0.));
    m_ikTree.InsertLeftChild(m_ikNodes[3], m_ikNodes[4]);
    
    m_ikNodes[5] = new Node(VectorR3(0, -0.164751, 0.996802), unity, 0.08, JOINT, minTheta, maxTheta, ofDegToRad(0.));
    m_ikTree.InsertLeftChild(m_ikNodes[4], m_ikNodes[5]);
    
    m_ikNodes[6] = new Node(VectorR3(0, -0.164751, 0.996802), zero, 0.08, EFFECTOR, minTheta, maxTheta, ofDegToRad(0.));
    m_ikTree.InsertLeftChild(m_ikNodes[5], m_ikNodes[6]);
    
//    //const VectorR3& unitx = VectorR3::UnitX;
//    const VectorR3& unity = VectorR3::UnitY;
//    const VectorR3& unitz = VectorR3::UnitZ;
//    const VectorR3 unit1(sqrt(14.0) / 8.0, 1.0 / 8.0, 7.0 / 8.0);
//    const VectorR3& zero = VectorR3::Zero;
//    
//    float minTheta = -4 * PI;
//    float maxTheta = 4 * PI;
//    
//    m_ikNodes.resize(8);//7DOF+additional endeffector
//    
//    m_ikNodes[0] = new Node(VectorR3(0.100000, 0.000000, 0.087500), unitz, 0.08, JOINT, -1e30, 1e30, ofDegToRad(0.));
//    m_ikTree.InsertRoot(m_ikNodes[0]);
//    
//    m_ikNodes[1] = new Node(VectorR3(0.100000, -0.000000, 0.290000), unity, 0.08, JOINT, -0.5, 0.4, ofDegToRad(0.));
//    m_ikTree.InsertLeftChild(m_ikNodes[0], m_ikNodes[1]);
//    
//    m_ikNodes[2] = new Node(VectorR3(0.100000, -0.000000, 0.494500), unitz, 0.08, JOINT, minTheta, maxTheta, ofDegToRad(0.));
//    m_ikTree.InsertLeftChild(m_ikNodes[1], m_ikNodes[2]);
//    
//    m_ikNodes[3] = new Node(VectorR3(0.100000, 0.000000, 0.710000), -unity, 0.08, JOINT, minTheta, maxTheta, ofDegToRad(0.));
//    m_ikTree.InsertLeftChild(m_ikNodes[2], m_ikNodes[3]);
//    
//    m_ikNodes[4] = new Node(VectorR3(0.100000, 0.000000, 0.894500), unitz, 0.08, JOINT, minTheta, maxTheta, ofDegToRad(0.));
//    m_ikTree.InsertLeftChild(m_ikNodes[3], m_ikNodes[4]);
//    
//    m_ikNodes[5] = new Node(VectorR3(0.100000, 0.000000, 1.110000), unity, 0.08, JOINT, minTheta, maxTheta, ofDegToRad(0.));
//    m_ikTree.InsertLeftChild(m_ikNodes[4], m_ikNodes[5]);
//    
//    m_ikNodes[6] = new Node(VectorR3(0.100000, 0.000000, 1.191000), unitz, 0.08, JOINT, minTheta, maxTheta, ofDegToRad(0.));
//    m_ikTree.InsertLeftChild(m_ikNodes[5], m_ikNodes[6]);
//    
//    m_ikNodes[7] = new Node(VectorR3(0.100000, 0.000000, 1.20000), zero, 0.08, EFFECTOR);
//    m_ikTree.InsertLeftChild(m_ikNodes[6], m_ikNodes[7]);
    
}


//--------------------------------------------------------------
void ofApp::keyPressed(int key){
}

//--------------------------------------------------------------
void ofApp::keyReleased(int key){
    
}

//--------------------------------------------------------------
void ofApp::mouseMoved(int x, int y ){
    
}

//--------------------------------------------------------------
void ofApp::mouseDragged(int x, int y, int button){
    
}

//--------------------------------------------------------------
void ofApp::mousePressed(int x, int y, int button){
    
}

//--------------------------------------------------------------
void ofApp::mouseReleased(int x, int y, int button){
    
}

//--------------------------------------------------------------
void ofApp::mouseEntered(int x, int y){
    
}

//--------------------------------------------------------------
void ofApp::mouseExited(int x, int y){
    
}

//--------------------------------------------------------------
void ofApp::windowResized(int w, int h){
    
}

//--------------------------------------------------------------
void ofApp::gotMessage(ofMessage msg){
    
}

//--------------------------------------------------------------
void ofApp::dragEvent(ofDragInfo dragInfo){
    
}
//
