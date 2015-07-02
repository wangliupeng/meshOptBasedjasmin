//
// 文件名: main_MeshOpt.C
//

#include "JASMIN_config.h"

#include "tbox/JASMINManager.h"
#include "tbox/MPI.h"
#include "tbox/InputDatabase.h"
#include "tbox/InputManager.h"
#include "tbox/RestartManager.h"
#include "tbox/Timer.h"
#include "tbox/TimerManager.h"
#include "tbox/Utilities.h"
#include "tbox/MemoryUtilities.h"

#include "MultiblockPatchHierarchy.h"
#include "CartesianCoordinates.h"
#include "MultiblockDeformingGridGeometry.h"
#include "HierarchyTimeIntegrator.h"
#include "DomainDataWriter.h"

#include "DeformingGridInputUtilities2.h"

#include "MeshOptLevelIntegrator.h"
#include "MeshOpt.h"

using namespace JASMIN;


int main( int argc, char *argv[] )
 { 

  /*************************************************************************
   *                                                                       *
   * 一、初始化.                                                           *
   * 包括: JASMIN环境初始化、读入输入参数、                                *
   * 创建应用中需要的JASMIN类对象和自定义类对象、                          *
   * 初始化网格片层次结构和数据片数据等等.                                 *
   *                                                                       *
   *************************************************************************/

   // 初始化MPI和JASMIN环境.      
   tbox::MPI::init(&argc, &argv);
   tbox::JASMINManager::startup();

   {
   // 解析命令行参数. 获取计算参数的输入文件名, 
   string input_filename;
   if ( (argc != 2) && (argc != 4) ) {
        tbox::pout << "USAGE:  " << argv[0] << " <input filename> "
           << "<restart dir> <restore number> [options]\n"
           << "  options:\n"
           << "  none at this time"
           << endl;
        tbox::MPI::abort();
        return (-1);
   } else {
      input_filename = argv[1];   
   }

   // 创建并解析输入文件的计算参数到输入数据库, 称之为根数据库.
   tbox::Pointer<tbox::Database> input_db = new tbox::InputDatabase("input_db");
   tbox::InputManager::getManager()->parseInputFile(input_filename, input_db);


   // 创建计时器.
   tbox::TimerManager::createManager(input_db->getDatabase("TimerManager"));
   
   
  /**************************************************************************/

   // 创建JASMIN框架类的对象, 以及自定义类对象, 并初始化. 
   
   // 创建笛卡尔坐标系.
   tbox::Pointer<geom::CoordinateSystem<NDIM> > coordinate_system = 
       new  geom::CartesianCoordinates<NDIM>();

   // 创建变形结构网格几何.
   tbox::Pointer<geom::MultiblockDeformingGridGeometry<NDIM> > grid_geometry =
       new geom::MultiblockDeformingGridGeometry<NDIM>("MultiblockGridGeometry",
                             coordinate_system); 

   //网格生成辅助工具类
   tbox::Pointer< appu::DeformingGridInputUtilities2 > grid_tool =
          new appu::DeformingGridInputUtilities2(
                     "DeformingGridInputUtilities2",
                     input_db->getDatabase("DeformingGridInputUtilities2"),
                     grid_geometry);

   // 创建网格片层次结构. 
   tbox::Pointer<hier::MultiblockPatchHierarchy<NDIM> > patch_hierarchy = 
       new hier::MultiblockPatchHierarchy<NDIM>("MultiblockPatchHierarchy", 
                                                grid_geometry);

                                                   
   // 创建网格片上MeshOpt类对象
   MeshOpt* MeshOpt_model = new MeshOpt("MeshOpt", 
                                          input_db->getDatabase("MeshOpt"), 
                                          grid_tool,
                                          grid_geometry);    	
                                                                  
   //7.5 创建网格层时间积分算法.
   MeshOptLevelIntegrator *MeshOpt_level_integrator =
	 new MeshOptLevelIntegrator("MeshOptLevelIntegrator",
	                input_db->getDatabase("MeshOptLevelIntegrator"),
	                MeshOpt_model,
                        grid_geometry);  


   //7.6 网格片层次结构显式时间积分算法类.
   tbox::Pointer<algs::HierarchyTimeIntegrator<NDIM> > time_integrator =
	 new algs::HierarchyTimeIntegrator<NDIM>("HierarchyTimeIntegrator",
			input_db->getDatabase("HierarchyTimeIntegrator"),
			patch_hierarchy,
			MeshOpt_level_integrator);

  


 

 /**********************************************************************/

      // *** 初始化网格片层次结构和所有网格片上的数据. ***    
   if(!time_integrator->initializeHierarchy()) {
         TBOX_ERROR("\nHierarchy is not successfully initialized. "<<endl);
      }
   
   
  

  /*************************************************************************
   *                                                                       *
   *  二、主循环                                                           *
   *                                                                       *
   *************************************************************************/

   double loop_time = time_integrator->getIntegratorTime();
   double loop_time_end = time_integrator->getEndTime();

   while ( (loop_time < loop_time_end) && 
          time_integrator->stepsRemaining() ) {

      int iteration_num = time_integrator->getIntegratorStep() + 1;

      tbox::plog <<endl<<endl;
      tbox::pout << "++++++++++++++++++++++++++++++++++++++++++++" << endl;
      tbox::pout << "At begining of timestep # " << iteration_num - 1 << endl;
      tbox::pout << "Simulation time is " << loop_time << endl;
      tbox::pout << "++++++++++++++++++++++++++++++++++++++++++++\n\n" << endl;
      tbox::plog <<endl<<endl;

         // 将数值解推进一个时间步, 返回其时间步长. 
	 double dt_actual = time_integrator->advanceHierarchy(); 

	 loop_time += dt_actual;

      tbox::plog <<endl<<endl;
      tbox::pout << "\n\n++++++++++++++++++++++++++++++++++++++++++++" << endl;
      tbox::pout << "At end of timestep # " <<  iteration_num - 1 << endl;
	 tbox::pout << "Dt = " << dt_actual << ", Simulation time is " 
                    << loop_time << endl;
      tbox::pout << "++++++++++++++++++++++++++++++++++++++++++++\n\n" << endl;
      tbox::plog <<endl<<endl; 

   }

  /*************************************************************************
   *                                                                       *
   *  三、模拟结束. 打印结果，并释放资源.                                  *
   *                                                                     *
   *************************************************************************/

   // 输出计时器统计的时间结果.
   tbox::TimerManager::getManager()->print(tbox::plog);

   tbox::MemoryUtilities::printMemoryInfo(tbox::plog);
   tbox::MemoryUtilities::printMaxMemory(tbox::plog);

   // 释放标准指针.
   if (MeshOpt_model) delete MeshOpt_model;
   if (MeshOpt_level_integrator) delete MeshOpt_level_integrator;

   }

   // 注销JASMIN系统.
   tbox::JASMINManager::shutdown();
   tbox::MPI::finalize();

   return(0); 
}

