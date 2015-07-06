 //
// 文件名:  MeshOptLevelIntegrator.h
//

#ifndef included_MeshOptLevelIntegrator
#define included_MeshOptLevelIntegrator

#include "JASMIN_config.h"
#include "BasePatchHierarchy.h"
#include "BasePatchLevel.h"
#include "tbox/Pointer.h"

#include <set>
using namespace std;

#include "StandardComponentPatchStrategy.h"
#include "TimeIntegratorLevelStrategy.h"

#include "NumericalIntegratorComponent.h"
#include "DtIntegratorComponent.h"
#include "InitializeIntegratorComponent.h"
#include "SynchronizeIntegratorComponent.h"
#include "CopyIntegratorComponent.h"
#include "MemoryIntegratorComponent.h"
#include "ReductionIntegratorComponent.h"
#include "OuterdataOperationIntegratorComponent.h"

#include "MultiblockGridGeometry.h"

using namespace JASMIN;

#include "MeshOpt.h"

/**
 * @brief 该类实现网格层时间积分算法策略类.
 * 该类仅考虑单层网格.
 *
 * 该类实现如下抽象接口函数:
 *    -# initializeLevelIntegrator()
 *    -# initializeLevelData()
 *    -# getLevelDt()
 *    -# advanceLevel()
 *    -# acceptTimeDependentSolution()
 *
 *
 * @see algs::TimeIntegratorLevelStrategy, algs::IntegratorComponentManager
 *
 */

class MeshOptLevelIntegrator :
public tbox::Serializable,
public algs::TimeIntegratorLevelStrategy<NDIM>
{
public:
   /**
    * @brief 构造函数.
    */
   MeshOptLevelIntegrator(const string& object_name,
           tbox::Pointer<tbox::Database> input_db,
           MeshOpt* patch_strategy,
           tbox::Pointer< hier::MultiblockGridGeometry<NDIM> > grid_geometry);

   /**
    * @brief 析构函数.
    */
   virtual ~MeshOptLevelIntegrator();


   /**
    * @brief 初始化网格层时间积分算法.
    */
   void
   initializeLevelIntegrator(tbox::Pointer<algs::IntegratorComponentManager<NDIM> > manager);


   /**
    * @brief 初始化指定网格层的数据片.
    * 初值构件类 algs::InitializeIntegratorComponent 支撑该函数的实现.
    */
   void
   initializeLevelData(const tbox::Pointer< hier::BasePatchHierarchy<NDIM> > hierarchy,
                       const int    level_number,
                       const double init_data_time,
                       const bool   can_be_refined,
                       const bool   initial_time,
                       const tbox::Pointer< hier::BasePatchLevel<NDIM> > old_level =
                               tbox::Pointer< hier::BasePatchLevel<NDIM> >(NULL),
                       const bool allocate_data = true);

   /**
    * @brief 获取指定网格层的时间步长. 步长构件类 algs::DtIntegratorComponent
    * 支撑该函数的实现.
    */
   double
   getLevelDt(const tbox::Pointer< hier::BasePatchLevel<NDIM> > level,
              const double dt_time,
              const bool initial_time,
              const int  flag_last_dt,
              const double last_dt);

   /**
    * @brief 网格层积分一个时间步. 数值构件类 algs::NumericalIntegratorComponent
    * 支撑该函数的实现.
    *
    */
   int
   advanceLevel(const tbox::Pointer< hier::BasePatchLevel<NDIM> > level,
                const double current_time,
                const double predict_dt,
                const double max_dt,
                const double min_dt,
                const bool   first_step,
                const int    hierarchy_step_number,
                double&      actual_dt);

   /**
    * @brief 在指定的网格层上, 时间步积分完毕, 存储最新的获得数值解,
    * 更新网格层的状态到新的时刻, 完成下一个时间步积分的准备工作.
    * 复制构件 algs::CopyIntegratorComponent 支撑该函数的实现.
    *
    */
   void
   acceptTimeDependentSolution(
                  const tbox::Pointer< hier::BasePatchLevel<NDIM> > level,
                  const double new_time,
                  const bool deallocate_data);

  void putToDatabase(tbox::Pointer<tbox::Database> db);

private:
   void getFromInput(
    tbox::Pointer<tbox::Database> db);

   // 对象名称.
   string d_object_name;


   // 网格片时间积分算法对象.
   MeshOpt* d_patch_strategy;

   // 网格几何对象
   tbox::Pointer< hier::MultiblockGridGeometry<NDIM> > d_grid_geometry;

   // 初值构件.
   tbox::Pointer< algs::InitializeIntegratorComponent<NDIM> >
                                                       d_init_intc;

   // 步长构件.
   tbox::Pointer< algs::DtIntegratorComponent<NDIM> > d_dt_intc;

   // 数值构件: 输出网格, 扰动网格, 优化网格.
   tbox::Pointer< algs::NumericalIntegratorComponent<NDIM> >
                                                   d_write_mesh_intc,
                                                   d_distrub_mesh_intc,
                                                   d_optimize_mesh_intc;
    
    // 外表面操作积分构件.
   tbox::Pointer< algs::OuterdataOperationIntegratorComponent<NDIM> >
                                                     d_outer_data_intc;

   // 复制构件: 移动网格迭代.
   tbox::Pointer< algs::CopyIntegratorComponent<NDIM> >
                                                       d_reset_intc,
                                                       d_cycle_copy_intc;


   // 内存构件: 新值, 演算.
   tbox::Pointer< algs::MemoryIntegratorComponent<NDIM> >
                                                       d_new_intc,
                                                       d_scratch_intc;

};

#endif
