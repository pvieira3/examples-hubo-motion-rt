#include "walker.h"

#define HAVE_HUBO_ACH
#include <hubo-zmp.h>

#include "HuboKin.h"
#include <string>

ach_channel_t zmp_chan;



/*
const double nudgePGain = 0.04;
const double nudgeIGain = 0.2;
const double nudgeDGain = 0.0;
*/


// gain for ankle flat
const double k_tau = 0.0010; // Grey likes this value

// gain for moment reducing
const double k_pos = 0.0030; // Crane stance resistance to torque

// gain for IMU rotation
const double k_theta_x = 0;
const double k_theta_y = 0.5*M_PI/180.0; // IMU feedback gain

const double weight_thresh_N = -1e5; // TESTING k_tau 30; // weight to turn on/off integrators
const double nudge_max_norm = 0.05; // m
const double spin_max_angle = 30 * M_PI/180;
const double comp_max_angle = 30 * M_PI/180;

const double tau_dead_band = 1;


const double fzMin = 10;
const double fzMax = 50;

const double flatteningGain = 0.02;


//const double compPitchGain = 0.0015;
//const double shiftCompMultiplier = 1.0;
//const double craneCompMultiplier = 1.0;


// TODO: resurrect this later when roll ankle limits give us trouble
//
//const double craneHipPitchGain = 0.0*M_PI/180.0; // Not used
//const double craneHipRollGain  = 0.05*M_PI/180.0; // Used when ankle roll reaches joint limits
//const double craneAnkleRollThresh = 2.0*M_PI/180.0;

//const double craneHipSafetyMargin = 2.5*M_PI/180.0; // Buffer to prevent leg collisions during crane stance
//const double craneHipSafetyGain = 1.5;

//const double swingVelLimit = 0.05;

//const double rollAngleGain = 0.3*M_PI/180.0; // IMU feedback gain
//const double craneAngleMultiplier = 1.0;

//const double shiftGain = 0.0002; // Shift-to-side gain
//const double shiftBalanceGain = 0.0015; // Gain to maintain balance while shifting
//const double craneTransitionPercent = 98.0; // Percent of weight on swing foot before entering crane stance

//const double crouchGain = 0.30;  // Leg-lifting gain

//const double pcAccel = 0.4; // Nominal acceleration used for position control
//const double vcAccel = 1.5; // Nominal acceleration used for velocity control

//const double calKneeGain = 1.0;
//const double calWeightGain = 0.001; // rad/sec per N

double applyDeadband( double x ) {

  if (x > tau_dead_band) {
    return x - tau_dead_band;
  } else if (x < -tau_dead_band) {
    return x + tau_dead_band;
  } else {
    return 0;
  }

}


void getLegNudge( Hubo_Control &hubo, Vector6d q, 
		  zmp_traj_element_t elem, nudge_state_t &state, 
		  int side, double dt )
{
    Eigen::Vector3d vel; vel.setZero();

    double refAngleX=0, refAngleY=0;
    // TODO: Calculate desired IMU angle

    Eigen::Matrix3d R = ( Eigen::AngleAxisd(q(HY), Vector3d::UnitZ()) *
			  Eigen::AngleAxisd(q(HR), Vector3d::UnitX()) *
			  Eigen::AngleAxisd(q(HP), Vector3d::UnitY()) *
			  Eigen::AngleAxisd(q(KN), Vector3d::UnitY()) *
			  Eigen::AngleAxisd(q(AP), Vector3d::UnitY()) *
			  Eigen::AngleAxisd(q(AR), Vector3d::UnitX()) ).toRotationMatrix();

    /*
    std::cout << "left MX: " << hubo.getLeftFootMx() << ", "
	      << "desired: " << elem.torque[LEFT][0] << ", "
	      << "diff: " << (hubo.getLeftFootMx() - elem.torque[LEFT][0]) << "\n";

    std::cout << "left MY: " << hubo.getLeftFootMy() << ", "
	      << "desired: " << elem.torque[LEFT][1] << ", "
	      << "diff: " << (hubo.getLeftFootMy() - elem.torque[LEFT][1]) << "\n";

    std::cout << "right MX: " << hubo.getRightFootMx() << ", "
	      << "desired: " << elem.torque[RIGHT][0] << ", "
	      << "diff: " << (hubo.getRightFootMx() - elem.torque[RIGHT][0]) << "\n";

    std::cout << "right MY: " << hubo.getRightFootMy() << ", "
	      << "desired: " << elem.torque[RIGHT][1] << ", "
	      << "diff: " << (hubo.getRightFootMy() - elem.torque[RIGHT][1]) << "\n\n";
    */

    if( side==RIGHT )
    {
      vel(0) =  k_pos*applyDeadband( hubo.getRightFootMy() - elem.torque[RIGHT][1] );
      vel(1) = -k_pos*applyDeadband( hubo.getRightFootMx() - elem.torque[RIGHT][0] );


        state.nudge += R*vel*dt;

        state.spin(0) += dt*k_theta_x*(hubo.getAngleX()-refAngleX);//-state.imu_offset(0));
        state.spin(1) += dt*k_theta_y*(hubo.getAngleY()-refAngleY);//-state.imu_offset(0));

        state.rarerr += dt*k_tau*( hubo.getRightFootMx() );
        state.raperr += dt*k_tau*( hubo.getRightFootMy() );
    }
    else
    {
        vel(0) =  k_pos*applyDeadband( hubo.getLeftFootMy() - elem.torque[LEFT][1] );
        vel(1) = -k_pos*applyDeadband( hubo.getLeftFootMx() - elem.torque[LEFT][0] );

        state.nudge += R*vel*dt;

        state.spin(0) += dt*k_theta_x*(hubo.getAngleX()-refAngleX);//-state.imu_offset(0));
        state.spin(1) += dt*k_theta_y*(hubo.getAngleY()-refAngleY);//-state.imu_offset(0));

        state.larerr += dt*k_tau*( hubo.getLeftFootMx() );
        state.laperr += dt*k_tau*( hubo.getLeftFootMy() );
    }



}

static inline void clamp(double& x, double cap) {
  x = std::max(-cap, std::min(x, cap));
}

// input is a hubo_control
// qr & ql are planned joint angles - used for 
// state holds the nudge integrators
// elem is our zmp_traj_element_t
// dt is likely to be 1.0/TRAJ_FREQ_HZ
void integrateNudgeIK( Hubo_Control &hubo, Vector6d &qr, Vector6d& ql,
		       nudge_state_t& state, zmp_traj_element_t elem, double dt,
		       HK::HuboKin& hkin) {

  // now truncate stuff
  if (state.nudge.norm() > nudge_max_norm) {
    state.nudge *= nudge_max_norm / state.nudge.norm();
  }

  clamp(state.spin[0], spin_max_angle);
  clamp(state.spin[1], spin_max_angle);
  clamp(state.larerr,  comp_max_angle);
  clamp(state.laperr,  comp_max_angle);
  clamp(state.rarerr,  comp_max_angle);
  clamp(state.raperr,  comp_max_angle);


  Vector6d* qboth[2] = { &qr, &ql };

  for (int side=0; side<2; ++side) {


    // Matt being paranoid: don't run IK code if no nudging to be done.
    if (state.nudge.squaredNorm() || state.spin.squaredNorm()) {

      Vector6d& qs = *(qboth[side]);


      Eigen::Isometry3d FK;

      // right comes first in HuboKin cause it was written by right handed people?
      hkin.legFK(FK, qs, side);

      // FK maps foot frame to body
      // 

      const Vector6d qs_prev = qs;

      const Eigen::Vector3d p = -state.nudge;
    
      const Eigen::Quaterniond R = 
	Eigen::AngleAxisd(state.spin[0], Eigen::Vector3d::UnitX()) *
	Eigen::AngleAxisd(state.spin[1], Eigen::Vector3d::UnitY());

      /*
      std::cout << "nudge for leg " << side << " is " << p.transpose() << "\n";
      std::cout << "spin for leg " << side << " is " << state.spin.transpose() << "\n";
      */

      Eigen::Isometry3d Tdelta;
      Tdelta.setIdentity();
      Tdelta.rotate(R);
      Tdelta.pretranslate(p);
      
      /*
      std::cout << "Tdelta for leg " << side << "=\n" << Tdelta.matrix() << "\n";
      Eigen::Isometry3d desired = Tdelta * FK;
      */

      Eigen::Isometry3d desired = FK;

      hkin.legIK(qs, desired, qs_prev, side);

      //std::cout << "Diff for leg " << side << " is " << (qs_prev - qs).transpose() << "\n";


    }

  }

  qr(AR) += state.rarerr;
  qr(AP) += state.raperr;

  ql(AR) += state.larerr;
  ql(AP) += state.laperr;


}

/*
// input is a hubo_control
// qr & ql are planned joint angles - used for 
// state holds the nudge integrators
// elem is our zmp_traj_element_t
// dt is likely to be 1.0/TRAJ_FREQ_HZ
void integrateNudge( Hubo_Control &hubo, Vector6d &qr, Vector6d& ql,
                        nudge_state_t state, zmp_traj_element_t elem, double dt )
{
    Vector6d dt_qr, dt_ql;
    int num = 20;
    double ddt = dt/num;
    

    for(int i=0; i<num; i++)
    {
        hubo.hipVelocityIK( dt_qr, state.nudge, state.spin, qr );
        hubo.hipVelocityIK( dt_ql, state.nudge, state.spin, ql );

        if( elem.stance == DOUBLE_LEFT || DOUBLE_RIGHT )
        {
            if( qr(AR)+dt_qr(AR)*ddt <= hubo.getJointAngleMin(RAR)
                || ( qr(HR)+qr(AR)<0.0 && qr(AR)<0 && qr(HR)>0 ) )
            {
                Eigen::Vector3d ghat( -sin(qr(HY)), cos(qr(HY)), 0 );
                Eigen::Vector3d altNudge = state.nudge - (state.nudge.dot(ghat))*ghat;

                hubo.hipVelocityIK( dt_qr, altNudge, state.spin, qr );
                hubo.hipVelocityIK( dt_ql, altNudge, state.spin, ql );
            }
            
            if( ql(AR)+dt_ql(AR)*ddt >= hubo.getJointAngleMax(LAR)
                || ( ql(HR)+ql(AR)>0.0 && ql(AR)>0 && ql(HR)<0 ) ) 
            {
                Eigen::Vector3d ghat( -sin(ql(HY)), cos(ql(HY)), 0 );
                Eigen::Vector3d altNudge = state.nudge - (state.nudge.dot(ghat))*ghat;

                hubo.hipVelocityIK( dt_qr, altNudge, state.spin, qr );
                hubo.hipVelocityIK( dt_ql, altNudge, state.spin, ql );
            }

            qr += dt_qr*dt;
            ql += dt_ql*dt;

        }
        else if( elem.stance == SINGLE_RIGHT )
        {
            if( qr(AR)+dt_qr(AR)*ddt <= hubo.getJointAngleMin(RAR)
                || ( qr(HR)+qr(AR)<0.0 && qr(AR)<0 && qr(HR)>0 ) )
            {
                Eigen::Vector3d ghat( -sin(qr(HY)), cos(qr(HY)), 0 );
                Eigen::Vector3d altNudge = state.nudge - (state.nudge.dot(ghat))*ghat;
                Eigen::Vector3d altSpin  = state.spin + craneHipRollGain/k_pos*
                                            ( (state.nudge.dot(ghat))*(ghat.cross(Eigen::Vector3d(0,0,1))) );
                state.imu_offset += (altSpin - state.spin)*ddt; 
                
                hubo.hipVelocityIK( dt_qr, altNudge, altSpin, qr );
                hubo.hipVelocityIK( dt_ql, altNudge, altSpin, ql );
            }

            qr += dt_qr*dt;
            ql += dt_ql*dt;

        }
        else if( elem.stance == SINGLE_LEFT )
        {
            if( ql(AR)+dt_ql(AR)*ddt >= hubo.getJointAngleMax(LAR)
                || ( ql(HR)+ql(AR)>0.0 && ql(AR)>0 && ql(HR)<0 ) )
            {
                Eigen::Vector3d ghat( -sin(qr(HY)), cos(qr(HY)), 0 );
                Eigen::Vector3d altNudge = state.nudge - (state.nudge.dot(ghat))*ghat;
                Eigen::Vector3d altSpin  = state.spin + craneHipRollGain/k_pos*
                                            ( (state.nudge.dot(ghat))*(ghat.cross(Eigen::Vector3d(0,0,1))) ); 
                state.imu_offset += (altSpin - state.spin)*ddt; 
                
                hubo.hipVelocityIK( dt_qr, altNudge, altSpin, qr );
                hubo.hipVelocityIK( dt_ql, altNudge, altSpin, ql );
            }

            qr += dt_qr*dt;
            ql += dt_ql*dt;

        }
        
    }

    qr(AR) += state.rarerr*ddt;
    qr(AP) += state.raperr*ddt;
    ql(AR) += state.larerr*ddt;
    ql(AP) += state.laperr*ddt;


}
*/

void nudgeRefs( Hubo_Control &hubo, zmp_traj_element_t &elem, //Eigen::Vector3d &vprev,
		nudge_state_t &state, double dt,
		HK::HuboKin& hkin)
{

  Vector6d qr, ql;
  qr(HY) = elem.angles[RHY];
  qr(HR) = elem.angles[RHR];
  qr(HP) = elem.angles[RHP];
  qr(KN) = elem.angles[RKN];
  qr(AP) = elem.angles[RAP];
  qr(AR) = elem.angles[RAR];

  ql(HY) = elem.angles[LHY];
  ql(HR) = elem.angles[LHR];
  ql(HP) = elem.angles[LHP];
  ql(KN) = elem.angles[LKN];
  ql(AP) = elem.angles[LAP];
  ql(AR) = elem.angles[LAR];

  // TODO: smoothly phase in the integrator for each foot
  if( hubo.getRightFootFz()+hubo.getLeftFootFz() > weight_thresh_N) {
        
      nudge_state_t lstate=state, rstate=state;

      getLegNudge( hubo, qr, elem, rstate, RIGHT, dt );
      getLegNudge( hubo, ql, elem, lstate, LEFT,  dt );
        
      state.larerr = lstate.larerr;
      state.rarerr = rstate.rarerr;
      state.laperr = lstate.laperr;
      state.raperr = rstate.raperr;

      if( elem.forces[RIGHT][2]+elem.forces[LEFT][2] == 0 ) {

	elem.forces[RIGHT][2] = 1;
	elem.forces[LEFT][2]  = 1;
	fprintf(stderr, "Warning: predicted forces both 0!");
      }

      state.nudge = (elem.forces[RIGHT][2]*rstate.nudge + elem.forces[LEFT][2]*lstate.nudge)/
	(elem.forces[RIGHT][2]+elem.forces[LEFT][2]);
      state.spin  = (elem.forces[RIGHT][2]*rstate.spin  + elem.forces[LEFT][2]*lstate.spin )/
	(elem.forces[RIGHT][2]+elem.forces[LEFT][2]);


  }


  std::cout << "Nudge: " << state.nudge.transpose() << "\tSpin: " << state.spin.transpose() << "\n";
  std::cout << "Offsets: " 
	    << state.raperr << ", " 
	    << state.rarerr << ", " 
	    << state.laperr << ", " 
	    << state.larerr << "\n";

  integrateNudgeIK( hubo, qr, ql, state, elem, dt, hkin );

        
  elem.angles[RHY] = qr(HY);
  elem.angles[RHR] = qr(HR);
  elem.angles[RHP] = qr(HP);
  elem.angles[RKN] = qr(KN);
  elem.angles[RAP] = qr(AP);
  elem.angles[RAR] = qr(AR);

  elem.angles[LHY] = ql(HY);
  elem.angles[LHR] = ql(HR);
  elem.angles[LHP] = ql(HP);
  elem.angles[LKN] = ql(KN);
  elem.angles[LAP] = ql(AP);
  elem.angles[LAR] = ql(AR);

}




void flattenFoot( Hubo_Control &hubo, zmp_traj_element_t &elem, nudge_state_t &state, double dt )
{
    
//    std::cout << "RFz:" << hubo.getRightFootFz() << "RAP:" << state.raperr << "\tLFz:" << hubo.getLeftFootFz() << "\tLAP:" << state.laperr << std::endl;

    if( fzMin < hubo.getRightFootFz() && hubo.getRightFootFz() < fzMax )
    {
        std::cout << "Flattening Right Foot" << std::endl;
        state.rarerr += dt*flatteningGain*( hubo.getRightFootMx() );
        state.raperr += dt*flatteningGain*( hubo.getRightFootMy() );
    }

    if( fzMin < hubo.getLeftFootFz() && hubo.getLeftFootFz() < fzMax )
    {
        std::cout<< "Flattening Left Foot" << std::endl;
        state.larerr += dt*flatteningGain*( hubo.getLeftFootMx() );
        state.laperr += dt*flatteningGain*( hubo.getLeftFootMy() );
    }

    elem.angles[RAR] += state.rarerr;
    elem.angles[RAP] += state.raperr;
    elem.angles[LAR] += state.larerr;
    elem.angles[LAP] += state.laperr;

}

void balance( Hubo_Control &hubo )
{
    std::cout << "balancing ... \n";

    double compRollGain = 0.0015;
    double compPitchGain = 0.0015;
    double pitchAngleGain = 0.5*M_PI/180.0;
    double rollAngleGain = 0.3*M_PI/180;

    double leftP=0, leftR=0, rightP=0, rightR=0; 

    double dt=0;
    double shutdownTime=5; // seconds
    double absoluteTime=0; // absolute time
    double time = hubo.getTime(); // initial time (s)
    int i=0, imax=40;

    // set ankle joints to velocity control
    hubo.setVelocityControl(RAP);
    hubo.setVelocityControl(RAR);
    hubo.setVelocityControl(LAP);
    hubo.setVelocityControl(LAR);

    while(absoluteTime < shutdownTime)
    {
        hubo.update(true);
        dt = hubo.getTime() - time;
        time = hubo.getTime();

        if(dt > 0)
        {
            absoluteTime += dt;

            // compute ankle joint velocities based on IMU and ankle torques
            leftP = pitchAngleGain*hubo.getAngleY() + compPitchGain*hubo.getLeftFootMy();
            leftR = rollAngleGain*hubo.getAngleX() + compRollGain*hubo.getLeftFootMx();
            rightP = pitchAngleGain*hubo.getAngleY() + compPitchGain*hubo.getRightFootMy();
            rightR = rollAngleGain*hubo.getAngleX() + compRollGain*hubo.getRightFootMx();

            // Set joint velocities
            hubo.setJointVelocity( RAP, rightP );
            hubo.setJointVelocity( RAR, rightR );
            hubo.setJointVelocity( LAP, leftP );
            hubo.setJointVelocity( LAR, leftR );

            // send controls
            //hubo.sendControls();

             // print output every imax cycles
            if( i>=imax )
            {
                std::cout //<< "\033[2J"
                          << "RAP Vel: " << rightP
                          << "\nRAR Vel: " << rightR
                          << "\nLAP Vel: " << leftP
                          << "\nLAR Vel: " << leftR
                          << "\nRt Torques: " << hubo.getRightFootMx() << ", " << hubo.getRightFootMy()
                          << "\nLt Torques: " << hubo.getLeftFootMx() << ", " << hubo.getLeftFootMy()
                          << std::endl;
            }
            if(i>=imax) i=0; i++;
        }
    }
}

// This function should be balancing and handling external forces too.
bool isStableCheck(Hubo_Control &hubo, double imuVelXInit, double imuVelYInit)
{
    bool isStable;
    double stableTol = 1.0;
    double rotVelX = hubo.getRotVelX();
    double rotVelY = hubo.getRotVelY();
    if(fabs(rotVelX-imuVelXInit) < stableTol && fabs(rotVelY-imuVelYInit) < stableTol)
        isStable = true;
    else
        isStable = false;
    return isStable;
}

/**
* @function: validateOutputData(TrajVector& traj)
* @brief: validation of joint angle output trajectory data
* @return: void
*/
void validateTrajSwapIn(double prevAngles[HUBO_JOINT_COUNT], double newAngles[HUBO_JOINT_COUNT]) {
    const double dt = 1.0/TRAJ_FREQ_HZ;
    double maxJointVel=0;
    double jointVel;
    const double jointVelTol = 6.0; // radians/s
    for (int j=0; j<HUBO_JOINT_COUNT; j++) {  
      jointVel = (prevAngles[j] - newAngles[j])/dt;
      if (jointVel > jointVelTol) {
        std::cerr << "change in joint " << jointNames[j] << " is larger than " << jointVelTol << "(" << jointVel << ")\n";
      }
      if (jointVel > maxJointVel) maxJointVel = jointVel;
    }
    std::cerr << "maxJntVel: " << maxJointVel << std::endl;
}


int main(int argc, char **argv)
{
    // create objects
    Hubo_Control hubo;
    HK::HuboKin hkin;

    hkin.kc.leg_l1 = 0; // eliminate neck -> waist Z distance
    hkin.kc.leg_l3 = 0; // eliminate waist -> hip Z distance
    hkin.kc.leg_l6 = 0; // eliminate ankle -> foot Z distance

    // create nudge state
    nudge_state_t state;
    memset( &state, 0, sizeof(nudge_state_t) );

    // open zmp ach channel
    ach_status_t r = ach_open( &zmp_chan, HUBO_CHAN_ZMP_TRAJ_NAME, NULL );
    fprintf( stderr, "ach_open: %s (%d)\n", ach_result_to_string(r), (int)r );

    // get zmp traj from ach channel
    size_t fs;
    zmp_traj_t curTrajectory, nextTrajectory;
    bool nextTrajReady = false;
    bool useNextTraj = false;
    walkState_t walkState = STOP;
    walkTransition_t walkTransition = STAY_STILL;
    bool isStable;

    // get initial rotational velocities of IMU
    double imuVelXInit = hubo.getRotVelX();
    double imuVelYInit = hubo.getRotVelY();

    // initialize time variables
    double dt, time, stime;
    double qDotTolerance = 0.005; // tolerance for joint change 
    double jointAngles[HUBO_JOINT_COUNT]; // to store current joint angles

    while(!daemon_sig_quit)
    {
        std::cout << "Getting trajectory\n";
        // if we don't have a new trajectory or we're told to stop,
        // then keep checking for walk trajectory
        while(curTrajectory.count <= 0 || curTrajectory.walkState == STOP)
        {
            memset( &curTrajectory, 0, sizeof(curTrajectory) );
            ach_status r = ach_get( &zmp_chan, &curTrajectory, sizeof(curTrajectory), &fs, NULL, ACH_O_LAST );
            if(r != ACH_STALE_FRAMES)
            {
                fprintf(stdout, "ach_get_result: %s\n", ach_result_to_string(r));
                std::cout << "count = " << curTrajectory.count << std::endl;
            }
        }

        fprintf(stderr, "Traj got. Count: %d\n", (int)curTrajectory.count);
        //for(int i=0; i<trajectory.count; i++)
        //    fprintf(stdout, "%d: RHR %f\n", i, trajectory.traj[i].angles[RHR] );

        // once we get a trajectory, check if we need to stabelize first or if we
        // can just keep walking, and then execute walk trajectory if not stopped
        if( curTrajectory.walkTransition == SWITCH_WALK )
        {
            std::cout << "stabilizing\n";
            while(isStable == false)
            {
                isStable = isStableCheck(hubo, imuVelXInit, imuVelYInit);
                std::cout << "stable! (hopefully)\n";
            }
        }

        if( curTrajectory.walkTransition == SWITCH_WALK )
        {
            // execute current trajectory
            for(int t=0; t<curTrajectory.count-1; t++)
            {
                // ####################################
                // ##### GET INTO WALK POSITION #######
                // ####################################
                if( t == 0 && curTrajectory.walkTransition == SWITCH_WALK )
                {
                    std::cout << "getting into walk position\n";
                    // update and set initial joint positions, speeds and accelerations
                    hubo.update(true);
                    for(int i=0; i<HUBO_JOINT_COUNT; i++)
                    {
                        hubo.setJointAngle( i, curTrajectory.traj[0].angles[i] );
                        hubo.setJointNominalSpeed( i, 0.4 );
                        hubo.setJointNominalAcceleration( i, 0.4 );
                    }

                    // set nominal speeds and acceleration of knee joints
                    hubo.setJointNominalSpeed( RKN, 0.8 );
                    hubo.setJointNominalAcceleration( RKN, 0.8 );
                    hubo.setJointNominalSpeed( LKN, 0.8 );
                    hubo.setJointNominalAcceleration( LKN, 0.8 );

                    // set initial should roll joint angles
                    hubo.setJointAngle( RSR, curTrajectory.traj[0].angles[RSR] + hubo.getJointAngleMax(RSR) );
                    hubo.setJointAngle( LSR, curTrajectory.traj[0].angles[LSR] + hubo.getJointAngleMin(LSR) );

                    // send commands
                    //hubo.sendControls();

                    // initialize times
                    stime=hubo.getTime(); time=hubo.getTime();

                    // update until 3 seconds has passed
                    while( time - stime < 3 ) {
                      hubo.update(true);
                      time = hubo.getTime();
                    }
                }

                // ####################################
                // ##### ELSE CONTINUE WALKING  #######
                // ####################################
                else
                {
                    std::cout << "continue walk, tick " << t << std::endl;
                    // if we don't have a next trajectory available
                    if( nextTrajReady == false )
                    {
                        // check for next trajectory
                        memset( &nextTrajectory, 0, sizeof(nextTrajectory) );
                        ach_get( &zmp_chan, &nextTrajectory, sizeof(nextTrajectory), &fs, NULL, ACH_O_LAST );

                        // if we received a next trajectory
                        if( nextTrajectory.trajNumber > curTrajectory.trajNumber )
                        {
                            // check if its start tick is already passed, indicate that the
                            // next trajectory is ready but we shouldn't use it
                            if( t > nextTrajectory.startTick || nextTrajectory.walkTransition == SWITCH_WALK )
                            {
                                useNextTraj = false;
                                nextTrajReady = true;
                            }
                            else // if we can use it, indicate that we can use it
                            {
                                useNextTraj = true;
                                nextTrajReady = true;
                            }
                        }
                    }
                    // if the next trajectory starts on the current tick,
                    // then reset t to 0 and start executing it from here
                    if( useNextTraj == true && nextTrajectory.startTick == t )
                    {
                        std::cout << "swapping in next trajectory\n";
                        // FIXME store current angles in order to check change next iteration
                        memset(&jointAngles, 0, sizeof(jointAngles));
                        memcpy(&jointAngles, &curTrajectory.traj[t].angles, sizeof(jointAngles));

                        curTrajectory = nextTrajectory;
                        t = 0;
                        nextTrajReady = false;
                    }

                    // validate trajectory swap in by checking angles before swap-in
                    // and after swap-in
                    if(t == 0)
                    {
                        std::cout << "validating\n";
                        validateTrajSwapIn(curTrajectory.traj[0].angles, jointAngles);
                    }

                    // update state, delta time and time
                    hubo.update(true);
                    dt = hubo.getTime() - time;
                    time = hubo.getTime();

                    // flatten feet with compliance
//uncomment before running                    flattenFoot( hubo, curTrajectory.traj[t], state, dt );
                    //nudgeRefs( hubo, trajectory.traj[t], state, dt, hkin ); //vprev, verr, dt );

                    // set joint angles for current trajectory tick
                    for(int i=0; i<HUBO_JOINT_COUNT; i++)
                    {
                        hubo.setJointAngle( i, curTrajectory.traj[t].angles[i] );
                        hubo.setJointNominalSpeed( i,
                               (curTrajectory.traj[t].angles[i]-curTrajectory.traj[t-1].angles[i])*TRAJ_FREQ_HZ );
                        double accel = TRAJ_FREQ_HZ*TRAJ_FREQ_HZ*(
                                            curTrajectory.traj[t-1].angles[i]
                                        - 2*curTrajectory.traj[t].angles[i]
                                        +   curTrajectory.traj[t+1].angles[i] );
                        hubo.setJointNominalAcceleration( i, 10*accel );
                    }

                    // set shoulder roll joint angles
                    hubo.setJointAngle( RSR, curTrajectory.traj[t].angles[RSR] + hubo.getJointAngleMax(RSR) );
                    hubo.setJointAngle( LSR, curTrajectory.traj[t].angles[LSR] + hubo.getJointAngleMin(LSR) );

                    // set hip roll joint angles
                    hubo.setJointAngleMin( LHR, curTrajectory.traj[t].angles[RHR] );
                    hubo.setJointAngleMax( RHR, curTrajectory.traj[t].angles[LHR] );

                    // send commands
                    //hubo.sendControls();
                }
            }
        }
    }
    return 0;
}
