//
// 文件名:	MeshOpt.h
//
 
#ifndef included_MeshOpt
#define included_MeshOpt

#ifndef included_JASMIN_config
#include "JASMIN_config.h"
#endif

#include <set>
using namespace std;

#include "tbox/Pointer.h"
#include "tbox/Array.h"
#include "tbox/Database.h"
#include "tbox/PIO.h"
#include "tbox/Serializable.h"

#include "CellVariable.h"
#include "NodeVariable.h"
#include "CellData.h"
#include "NodeData.h"
#include "SideVariable.h"
#include "SideData.h"
#include "VariableContext.h"
#include "Patch.h"

#include "DomainDataWriter.h"

#if (NDIM == 2)
#include "DeformingGridInputUtilities2.h"
#endif

#include "MultiblockDeformingGridGeometry.h"
#include "StandardComponentPatchStrategy.h"
#include "IntegratorComponent.h"

#include "Mesquite.hpp"
#include "Mesquite_all_headers.hpp"

using namespace Mesquite2;
using namespace JASMIN;


/**
 * @brief 类 MeshOpt 在单个网格片, 提供数值计算子程序, 
 * 
 * 实现网格片时间积分算法基类 algs::StandardComponentPatchStrategy.
 *
 * 实现如下纯虚函数:
 *  - initializeComponent() : 支撑所有积分构件的初始化.
 *  - computeOnpatch() :      支撑数值构件, 完成通信和计算.
 *  - getPatchDt()   :        支撑步长构件, 计算时间步长.
 *  - initializePatchData() : 为初值构件完成初始化.
 *
 *
 * @see algs::IntegratorComponent
 */

class MeshOpt : 
   public tbox::Serializable,
   public algs::StandardComponentPatchStrategy<NDIM>
{
public:

   MeshOpt(const string& object_name,
         tbox::Pointer<tbox::Database> input_db,
         tbox::Pointer< appu::DeformingGridInputUtilities2 > grid_tool,
         tbox::Pointer< geom::MultiblockDeformingGridGeometry<NDIM> > grid_geom);  

    /**
     * @brief 析构函数.
     */
   ~MeshOpt();

   
   /*!
    * @brief 初始化指定的积分构件,
    *        注册待填充的数据片或待调度内存空间的数据片到积分构件.
    */
     void initializeComponent(
            algs::IntegratorComponent<NDIM>* intc) const;

   /*!
    * @brief 当前积分构件为步长构件. 该函数计算时间步长.
    *
    */
    double getPatchDt(hier::Patch<NDIM>& patch,
                             const double  time,
                             const bool    initial_time,
                             const int     flag_last_dt,
                             const double  last_dt,
                             const string& intc_name);

   /*!
    * @brief 当前积分构件为数值构件. 该函数完成数值计算.
    *
    */
    void computeOnPatch(hier::Patch<NDIM>& patch,
                             const double  time,
                             const double  dt,
                             const bool    initial_time,
                             const string& intc_name);


   /**
    * @brief 当前积分构件为初值构件. 该函数初始化数据片.
    *
    */
    void initializePatchData(
      hier::Patch<NDIM>& patch,
      const double  time,
      const bool    initial_time,
      const string& intc_name);


   void setPhysicalBoundaryConditions(
                  hier::Patch<NDIM>& patch,
                  const double fill_time,
                  const hier::IntVector<NDIM>& ghost_width_to_fill,
                  const string& intc_name);


  void putToDatabase(tbox::Pointer<tbox::Database> db);
  

private:
  // 以下为私有成员函数.

   void getFromInput(
                tbox::Pointer<tbox::Database> db);

   // 创建所有变量和上下文, 注册变量上下文配对到变量数据库.
   void registerModelVariables();

   // 标定结点类型
   void setNodeInfo(hier::Patch<NDIM>& patch);

   //创建用于优化的网格数据结构
    MeshImpl * createLocalMesh(hier::Patch<NDIM> & patch);

  // 扰动网格
   void disturbMesh(hier::Patch<NDIM>& patch,
                    const double  time,
                    const double  dt,
                    const bool    initial_time);


   // 提取patch网格,存储到VTK文件
   void writeToVTK(hier::Patch<NDIM>& patch,
                   const double  time,
                   const double  dt,
                   const bool    initial_time);
 


   // 以下为私有数据成员.

    // 对象名.
   string d_object_name;     

   // 变形网格辅助输入工具.
   tbox::Pointer< appu::DeformingGridInputUtilities2 > d_grid_tool;

   //网格几何
   tbox::Pointer<geom::MultiblockDeformingGridGeometry<NDIM> > d_grid_geometry;



   // 存储坐标
   tbox::Pointer< pdat::NodeVariable<NDIM,double> > d_coords;

   // 固定点
   tbox::Pointer<pdat::NodeVariable<NDIM,bool> > d_fixed;
  

   // 上下文: 当前值, 新值, 演算值.
   tbox::Pointer< hier::VariableContext > d_current;
   tbox::Pointer< hier::VariableContext > d_new;
   tbox::Pointer< hier::VariableContext > d_scratch;
   tbox::Pointer< hier::VariableContext > d_move;

   // 数据片索引号.
   int d_coords_current_id,   // < d_coords, d_current >: 存储结点坐标的当前值.
       d_coords_new_id,       // < d_coords, d_new >    : 存储结点坐标的最新值.
       d_coords_scratch_id;   // < d_coords, d_scratch >: 存储结点坐标的演算值.

   // 数据片的两种影像区宽度.
   hier::IntVector<NDIM> d_zeroghosts, d_oneghosts;

    // 是否扰动
    bool d_flag;

    // 结点信息
    tbox::Array<bool> d_fixed_info;
                  
};

#endif
