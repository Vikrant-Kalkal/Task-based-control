#include <visp/vpFeaturePoint.h>
#include <ecn_baxter_vs/baxter_arm.h>
#include <visp/vpSubMatrix.h>
#include <visp/vpSubColVector.h>
#include <ecn_common/visp_utils.h>

using namespace std;

int main(int argc, char** argv)
{
    //BaxterArm arm(argc, argv);    // defaults to simulation with right arm
    BaxterArm arm(argc, argv, false, "left");   // real robot with left arm


    vpColVector q = arm.init();

    vpColVector qmin = arm.jointMin(),
                qmax = arm.jointMax();

    // define a simple 2D point feature and its desired value
    vpFeaturePoint p,pd;
    pd.set_xyZ(0,0,1);

    // the error
    vpColVector e(3);
    double x, y, area;
    // desired area
    const double area_d = arm.area_d();

    // loop variables
    vpColVector qdot(7), e_wave(10), q_star(7),q_actpostive(7), q_actnegative(7) ;
    vpMatrix L(3, 6), Js(3,7), H(10,10), J_wave(10, 7);
	q_star = arm.init();
	
    while(arm.ok())
        {
            cout << "-------------" << endl;

            // get point features
            x = arm.x();
            y = arm.y();
            area = arm.area();
            p.set_xyZ(x,y, 1);
            std::cout << "x: " << x << ", y: " << y << ", area: " << area << '\n';

            // update error vector e
            e[0] = x;
            e[1] = y;
            e[2] = area - area_d;

            // update interaction matrix L
			ecn::putAt(L, p.interaction(), 0, 0);
			L[2][0] = 0;
			L[2][1] = 0;
			L[2][2] = area*2;
			L[2][3] = 3*area*y;
			L[2][4] = -3*area*x;
			L[2][5] = 0;
            // compute feature Jacobian from L and cameraJacobian
			Js = L * arm.cameraJacobian();
			qdot = -arm.lambda()*Js.pseudoInverse()*e;
            // build H matrix (2nd section) using arm.rho()
            q = arm.jointPosition();
            q_star = (qmin + qmax)/2;
			
			double ro = arm.rho();	
			vpMatrix I3;
			I3.eye(3);
			vpMatrix I4;
			I4.eye(4);
			ecn::putAt(H, I3, 0, 0);
			ecn::putAt(H, I4, 6, 6);
			for(int i = 0; i<3; i++) {
				if (q[i] > q_star[i]) {
					q_actpostive[i] = qmax[i] - ro*(qmax[i] - qmin[i]);
					H[3+i][3+i] = ecn::weight(q[i], q_actpostive[i], qmax[i]);
				}else{
					q_actnegative[i] = qmin[i] + ro*(qmax[i] - qmin[i]);
					H[3+i][3+i] = ecn::weight(q[i], q_actnegative[i], qmin[i]);
				}		
			}
			
			ecn::putAt(e_wave, e, 0);
			ecn::putAt(e_wave, q - q_star, 3);
			
			vpMatrix I7;
			I7.eye(7);
			ecn::putAt(J_wave, Js, 0, 0);
			ecn::putAt(J_wave, I7, 3, 0);
			qdot = -arm.lambda()*(H*J_wave).pseudoInverse()*H*e_wave;
            // send this command to the robot
            arm.setJointVelocity(qdot);

            // display current joint positions and VS error
            arm.plot(e);
    }
}
