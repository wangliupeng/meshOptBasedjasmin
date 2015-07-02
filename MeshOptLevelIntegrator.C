 //
// 文件名:  MeshOptLevelIntegrator.C


#include <stdlib.h>
#include <fstream>
using namespace std;

#include "MeshOptLevelIntegrator.h" 
#include "MeshOpt.h" 

#include "MultiblockPatchLevel.h" 
#include "tbox/MPI.h" 

#ifdef DEBUG_CHECK_ASSERTIONS
#include <assert.h>
#endif

#include "tbox/IEEE.h"
#include "tbox/PIO.h"


using namespace JASMIN;

/****************************************************************
* 构造函数.
*****************************************************************/
MeshOptLevelIntegrator::MeshOptLevelIntegrator(
        const string& object_name,
        tbox::Pointer<tbox::Database> input_db,
        MeshOpt* patch_strategy,
        tbox::Pointer< hier::MultiblockGridGeometry<NDIM> > grid_geometry)
:
   d_object_name(object_name),
   d_patch_strategy(patch_strategy),
   d_grid_geometry(grid_geometry)
{
   getFromInput(input_db);
}

/****************************************************************
* 析构函数.
*****************************************************************/
MeshOptLevelIntegrator::~MeshOptLevelIntegrator() {}

/****************************************************************
*  创建并初始化所有积分构件.
*****************************************************************/
void MeshOptLevelIntegrator::initializeLevelIntegrator(
        tbox::Pointer<algs::IntegratorComponentManager<NDIM> > manager)
	{

	  // 初值构件 : 初始化新创建的网格层.
	  d_init_intc = new algs::InitializeIntegratorComponent<NDIM>("INIT",
                                                              d_patch_strategy,
                                                              true, manager);
	 
	  // 步长构件 : 计算时间步长.
	  d_dt_intc   = new algs::DtIntegratorComponent<NDIM>("TIME_STEP_SIZE",
                                                              d_patch_strategy,
                                                              manager);
	
	 // 数值构件: 写网格.
         d_mesh_intc = new algs::NumericalIntegratorComponent<NDIM>("WRITE_MESH",
                                                                 d_patch_strategy,
                                                                 manager); 
	 // 外表面操作积分构件.
	  d_outer_data_intc = new algs::OuterdataOperationIntegratorComponent<NDIM>("OUTER_DATA",
                                                              d_patch_strategy,
                                                              manager,
							      "MAX");
	 // 复制构件:  同步前复制数据片.
	  d_reset_intc   = new algs::CopyIntegratorComponent<NDIM>("RESET_SOLUTION",
                                                              d_patch_strategy,
                                                              manager);
	  d_cycle_copy_intc = new algs::CopyIntegratorComponent<NDIM>("COPY_PREMOVING",
                                                              d_patch_strategy,
                                                              manager);

	  // 内存构件: 当前值数据片, 新值数据片, 演算数据片.
	  d_new_intc     = new algs::MemoryIntegratorComponent<NDIM>("NEW_ALLOC_PATCH_DATA",
                                                              d_patch_strategy,
                                                              manager);
	  d_scratch_intc = new algs::MemoryIntegratorComponent<NDIM>("SCRATCH_ALLOC_PATCH_DATA",
                                                              d_patch_strategy,
                                                              manager);
	 

}

/****************************************************************
* 初始化网格层.
*****************************************************************/
void MeshOptLevelIntegrator::initializeLevelData(
      const tbox::Pointer< hier::BasePatchHierarchy<NDIM> > hierarchy, 
      const int    level_number,
      const double init_data_time,
      const bool   can_be_refined,
      const bool   initial_time,
      const tbox::Pointer< hier::BasePatchLevel<NDIM> > old_level,
      const bool allocate_data)
{
     NULL_USE(can_be_refined);

     // 初始化网格层.
     d_init_intc->initializeLevelData(hierarchy,
                                      level_number,
                                      init_data_time,
                                      initial_time,
                                      old_level,
                                      false);

}

/****************************************************************
* 计算时间步长.
*****************************************************************/
double MeshOptLevelIntegrator::getLevelDt(
              const tbox::Pointer< hier::BasePatchLevel<NDIM> > level, 
              const double dt_time, 
              const bool initial_time,
              const int  flag_last_dt,
              const double last_dt)
{
     double dt = d_dt_intc->getLevelDt(level,
                                       dt_time,
                                       initial_time,
                                       flag_last_dt,
                                       last_dt,
                                       false);

     return(dt);

}

/****************************************************************
* 积分一个时间步长.
*****************************************************************/
int MeshOptLevelIntegrator::advanceLevel(
                const tbox::Pointer< hier::BasePatchLevel<NDIM> > level, 
                const double current_time,
                const double predict_dt,
                const double max_dt,
                const double min_dt,
                const bool   first_step,
                const int    hierarchy_step_number,
                double&      actual_dt)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(predict_dt>=min_dt && predict_dt<=max_dt);
#endif

    NULL_USE(hierarchy_step_number);

    // 显式时间离散格式, 取时间步长等于预测时间步长.
    actual_dt = predict_dt;

    // 新值数据片调度内存空间.
    d_new_intc->allocatePatchData(level,current_time);
    d_scratch_intc->allocatePatchData(level,current_time);
    
    // 复制数据: current->new.
    d_cycle_copy_intc-> copyPatchData(level,
                                      current_time);

    //输出网格
    d_mesh_intc->computing(level,
                           current_time,
                           actual_dt); 

     // 外表面同步
    d_outer_data_intc->operate(level);        
                  
    // 设置新值数据片的时刻.
    d_new_intc->setTime(level,current_time+actual_dt);

    return(1);

}

/****************************************************************
* 接收数值解.
*****************************************************************/
void MeshOptLevelIntegrator::acceptTimeDependentSolution(
                const tbox::Pointer< hier::BasePatchLevel<NDIM> > level, 
                const double new_time,
                const bool deallocate_data)
{
    // 将数值解从新值数据片复制到当前值数据片.
    d_reset_intc->copyPatchData(level,
                                new_time);

    // 释放新值数据片.
    if(deallocate_data) d_new_intc->deallocatePatchData(level);

}

/****************************************************************
* 从输入文件中读取参数值.
*****************************************************************/
void MeshOptLevelIntegrator::getFromInput(
    tbox::Pointer<tbox::Database> db)
{

}

/****************************************************************
* 输出数据成员到重启动数据库. 
*****************************************************************/
void MeshOptLevelIntegrator::putToDatabase(
         tbox::Pointer<tbox::Database> db)
{

}


