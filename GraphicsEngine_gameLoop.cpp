#include "GraphicsEngine.h"

// ----------------------------------------------------------------------------
void GraphicsEngine::gameLoop(void)
{
	int iSimulationStatus = 0;
	while(e_GameState != GameState::EXIT)
	{
		drawGame();

		double dTimeIncrement_Request = 1.0e2*mpm_PhysicsEngine->getTime_Increment();
		if(iSimulationStatus == 0)
		{
//			iSimulationStatus = mpm_PhysicsEngine->runSimulation_Classic_SinglePass(dTimeIncrement_Request);
			iSimulationStatus = mpm_PhysicsEngine->runSimulation_Classic_SinglePass_MP(dTimeIncrement_Request);
//			iSimulationStatus = mpm_PhysicsEngine->runSimulation_Classic_SinglePass_MP_Contact(dTimeIncrement_Request);
		}

		processInput();
		this->f_Time = mpm_PhysicsEngine->getTime_Current();

		// save snapshot
		if(f_Time - f_TimeSnapshot_LastSave > f_TimeSnapshot_Interval)
		{
			drawGame();
			drawGame();

			f_TimeSnapshot_LastSave = f_Time;
			std::string strFileName = _STR_SNAPFILE;//"./bmp/Snapshot_";
			strFileName += Script(i_TimeSnapshotCycle);
			strFileName += ".bmp";
			this->saveScreenshot(0, 0, i_ScreenWidth, i_ScreenHeight, strFileName.c_str());
			i_TimeSnapshotCycle += 1;
		}
	}
}
// ----------------------------------------------------------------------------
