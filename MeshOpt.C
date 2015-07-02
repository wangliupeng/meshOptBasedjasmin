//
// 文件名:  MeshOpt.C
//

#ifndef included_MeshOpt_C
#define included_MeshOpt_C

#include "MeshOpt.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#include <assert.h>
#endif



#include "BlockPatchGeometry.h"
#include <iomanip>
#include "VariableDatabase.h"
#include "tbox/RestartManager.h"
#include "tbox/TimerManager.h"
#include "tbox/Timer.h"
#include "tbox/Utilities.h"
#include "tbox/IEEE.h"



MeshOpt::MeshOpt(const string& object_name,
         tbox::Pointer<tbox::Database> input_db,
         tbox::Pointer< appu::DeformingGridInputUtilities2 > grid_tool,
         tbox::Pointer<geom::MultiblockDeformingGridGeometry<NDIM> > grid_geom)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!object_name.empty());
   assert(!input_db.isNull());
   assert(!grid_geom.isNull());
   assert(!grid_tool.isNull());
#endif

   d_grid_geometry = grid_geom;
   d_grid_tool = grid_tool;

   getFromInput(input_db);

   d_flag = false;


   // 为影像区的宽度赋值.
   d_zeroghosts = hier::IntVector<NDIM>(0);
   d_oneghosts  = hier::IntVector<NDIM>(1);

   // 创建所有变量及数据片索引号, 注册可视化数据片.
   registerModelVariables();

}

/*
*************************************************************************
*    MeshOpt 类的空析构函数.	 	                                *
*************************************************************************
*/

MeshOpt::~MeshOpt() 
{ 

} 

/********************************************************************
* 创建变量和上下文, 注册变量上下文配对到变量数据库, 获得数据片索引号.
*********************************************************************/

void MeshOpt::registerModelVariables()
{
   // 变量数据库.
   hier::VariableDatabase<NDIM>* variable_db =
                  hier::VariableDatabase<NDIM>::getDatabase();

   // 定义变量. 
   d_coords   = new pdat::NodeVariable<NDIM,double>("coordinates", NDIM);
   
   // 当前值上下文, 新值上下文, 演算上下文, 可视化上下文.
   d_current = hier::VariableDatabase<NDIM>::getDatabase()->getContext("CURRENT");
   d_new     = hier::VariableDatabase<NDIM>::getDatabase()->getContext("NEW");
   d_scratch = hier::VariableDatabase<NDIM>::getDatabase()->getContext("SCRATCH");

   // 存储当前时刻的变量值: (coords,current), 影像区宽度为0.
   d_coords_current_id = 
         variable_db->registerVariableAndContext(d_coords,
                                                 d_current,
                                                 d_zeroghosts);

   // 存储新时刻的变量值: (coords,new), 影像区宽度为0.
   d_coords_new_id = 
         variable_db->registerVariableAndContext(d_coords,
                                                 d_new,
                                                 d_zeroghosts);

   // 存储变量的演算值: (coords,scratch), 影像区宽度等于格式宽度.
   d_coords_scratch_id = 
         variable_db->registerVariableAndContext(d_coords,
                                                 d_scratch,
                                                 d_oneghosts);

}


/********************************************************************
*  初始化积分构件.
*********************************************************************/
void MeshOpt::initializeComponent(algs::IntegratorComponent<NDIM>* intc) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(intc!=NULL);
#endif

   const string intc_name = intc->getName();

   if(intc_name=="INIT") {  // 初值构件 : 为网格坐标赋初值.
        intc->registerInitPatchData(d_coords_current_id);
   }else if(intc_name=="TIME_STEP_SIZE") { //步长构件: 求时间步长.

   }else if(intc_name=="COPY_PREMOVING") { //复制构件: 当前值赋值给新值.
        intc->registerCopyPatchData(d_coords_new_id,
                                d_coords_current_id);
   }else if(intc_name=="WRITE_MESH") { // 数值构件: 输出网格
	intc->registerRefinePatchData(d_coords_scratch_id,
                                      d_coords_new_id);
   }else if(intc_name=="DISTRUB_MESH") { // 数值构件: 输出网格
	intc->registerRefinePatchData(d_coords_scratch_id,
                                      d_coords_new_id);
   }else if(intc_name == "OUTER_DATA"){ // 外表面构件
          intc->registerPatchData (d_coords_new_id);
  }else if(intc_name=="RESET_SOLUTION") { // 复制构件: 接收数值解.
        intc->registerCopyPatchData(d_coords_current_id,
                                    d_coords_new_id);
        
   }else if(intc_name=="NEW_ALLOC_PATCH_DATA"){ // 内存构件: 为新值数据片调度内存空间.
        intc->registerPatchData(d_coords_new_id);
      
   }else if(intc_name=="SCRATCH_ALLOC_PATCH_DATA"){// 内存构件: 为演算数据片调度内存空间.
        intc->registerPatchData(d_coords_scratch_id);
      
   }else{
        TBOX_ERROR("\n::initializeComponent() : component "
                   << intc_name <<" is not matched. "<<endl);
   }    

}

/*
*************************************************************************
* 初值构件: 网格结点坐标
*************************************************************************
*/
void MeshOpt::initializePatchData( hier::Patch<NDIM>& patch,
                                 const double  time,
                                 const bool    initial_time,
                                 const string& intc_name)
{
   #ifdef DEBUG_CHECK_ASSERTIONS
   assert(intc_name=="INIT");
   #endif

   (void) time;

   if (initial_time) {
      tbox::Pointer< pdat::NodeData<NDIM,double> > coords_current = 
                               patch.getPatchData(d_coords_current_id);

#ifdef DEBUG_CHECK_ASSERTIONS
      assert(!coords_current.isNull());     
#endif

#ifdef DEBUG_CHECK_ASSERTIONS
      assert(ghost_cells == d_zeroghosts);
#endif

      // 生成初始网格.
      if(!d_grid_tool.isNull()) {
         d_grid_tool->generateDeformingMeshForDomain(patch,
                                                     d_coords_current_id);
      }

  }

}

/*
*************************************************************************
* 步长构件: 计算并返回网格片上的稳定时间步长.  
*************************************************************************
*/
double MeshOpt::getPatchDt(hier::Patch<NDIM>& patch,
                          const double  time,
                          const bool    initial_time,
                          const int     flag_last_dt,
                          const double  last_dt,
                          const string& intc_name)
{

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(intc_name=="TIME_STEP_SIZE");
#endif

   NULL_USE(flag_last_dt);
   NULL_USE(last_dt);

   tbox::Pointer< pdat::NodeData<NDIM,double> > coords_current = 
                               patch.getPatchData(d_coords_current_id);

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!coords_current.isNull());
  
#endif

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(ghost_cells == d_zeroghosts);
#endif

   
   double stabdt = 0.1;
  

   return stabdt;

}

/********************************************************************
*  实现数值构件.
*********************************************************************/
void MeshOpt::disturbMesh(hier::Patch<NDIM>& patch,
                    const double  time,
                    const double  dt,
                    const bool    initial_time)
{
    NULL_USE(dt);
    NULL_USE(initial_time);

    tbox::Pointer< pdat::NodeData<NDIM,double> > coords_current
                    = patch.getPatchData(d_coords_current_id);

    int count = -1;
    double dist = 0.01;
    for(pdat::NodeIterator<NDIM> ic((*coords_current).getBox()); ic; ic++)
    {
	dist *= -1;
        (*coords_current)(ic(),0) += dist;
        (*coords_current)(ic(),1) -= dist;
        ++count;
    }
    d_flag = true;

}


void MeshOpt::computeOnPatch(hier::Patch<NDIM>& patch,
                          const double  time,
                          const double  dt,
                          const bool    initial_time,
                          const string& intc_name)
{

   #ifdef DEBUG_CHECK_ASSERTIONS
   assert(intc_name=="WRITE_MESH"|| "DISTRUB_MESH");
   #endif

  if(intc_name =="WRITE_MESH") {      
       writeToVTK(patch,time,dt,initial_time);
    } else if(intc_name =="DISTRUB_MESH") {      
       disturbMesh(patch,time,dt,initial_time);
    }

}


/*
 ********************************************************************
 * 设置物理边界条件.                                                *
 ********************************************************************
 */ 
void MeshOpt::setPhysicalBoundaryConditions(
                  hier::Patch<NDIM>& patch,
                  const double fill_time,
                  const hier::IntVector<NDIM>& ghost_width_to_fill,
                  const string& intc_name)
{

}

void MeshOpt::writeToVTK(hier::Patch<NDIM>& patch,
                const double  time,
                const double  dt,
                const bool    initial_time)
{
    NULL_USE(dt);
    NULL_USE(initial_time);

   
         

    tbox::Pointer< pdat::NodeData<NDIM,double> > coords_current
                    = patch.getPatchData(d_coords_current_id);

    const hier::Index<NDIM> ifirst=patch.getBox().lower();
    const hier::Index<NDIM> ilast =patch.getBox().upper();

    int cols = ilast(0) - ifirst(0)+2;
    int rows = ilast(1) - ifirst(1)+2;

    int num_of_elems = (cols - 1)*(rows - 1);
    int num_of_nodes = cols*rows;
    tbox::Array<int> elem(4*num_of_elems,true);
    tbox::Array<double> node(3*num_of_nodes, true);

    int count = -1;
    for(int row = 0; row < rows-1; row++)
    {
        for(int col = 0; col < cols-1; col++)
        {
            int idx1 = row*cols+col;
            int idx2 = idx1 + 1;
            int idx3 = idx2 + cols;
            int idx4 = idx3 - 1;

            elem[++count] = idx1;
            elem[++count] = idx2;
            elem[++count] = idx3;
            elem[++count] = idx4;
        }

    }

    count = -1;
    for(pdat::NodeIterator<NDIM> ic((*coords_current).getBox()); ic; ic++)
    {
        node[++count]= (*coords_current)(ic(),0);
        node[++count]= (*coords_current)(ic(),1);
        node[++count]= 0.0;
    }

    //!todo test coords_current
    //!write to vtk

    const tbox::Pointer< hier::BlockPatchGeometry<NDIM> > pgeom =
                                            patch.getPatchGeometry();

    int block_index = pgeom->getBlockNumber();

    int patch_index = patch.getPatchNumber();

    std::stringstream bi, pi;
    bi << block_index;
    pi << patch_index;

    std::string file_name;
    if(d_flag)
	file_name = "1_block_ " + bi.str()+ "_patch_" +  pi.str()  + ".vtk";
    else 
        file_name = "0_block_ " + bi.str()+ "_patch_" +  pi.str()  + ".vtk";

    std::ofstream os(file_name.c_str());

    os<<"# vtk DataFile Version 3.0"<<"\n";
    os<<"CGAL mesh"<<"\n";
    os<<"ASCII"<<"\n";
    os<<"DATASET UNSTRUCTURED_GRID"<<"\n";

    //! 输出点的坐标
    os<<"POINTS " << num_of_nodes << " double" << "\n";
    os<< std::scientific << std::setprecision(12);
    for(int k = 0; k < num_of_nodes; k++)
    {
        os << node[3*k] << " "
           << node[3*k+1] << " "
           << node[3*k+2] << "\n";
    }

    //! 输出单元
    os<<"CELLS " << num_of_elems <<" "<< 5*num_of_elems <<"\n";
    for(int k = 0; k < num_of_elems; k++)
    {

        os<< "4 "
        << elem[4*k] <<" "
        << elem[4*k+1] <<" "
        << elem[4*k+2]<<" "
        << elem[4*k+3]<<"\n";

    }

    //! 输出单元类型
    os<<"CELL_TYPES "<< num_of_elems <<"\n";
    for(int i = 0; i < num_of_elems; i++)
    {
        os<< 9 <<"\n";
    }

    //! 输出节点属性
    os<<"POINT_DATA "<< num_of_nodes <<"\n";
    os<<"SCALARS fixed float"<<"\n"
     <<"LOOKUP_TABLE default"<<"\n";
    for(int k = 0; k < num_of_nodes; k++)
    {
        os << 1 << "\n";
    }

    //! 输出cell属性，属于block的编号
    os<<"CELL_DATA "<< num_of_elems <<"\n";
    os<<"SCALARS fixed float"<<"\n"
     <<"LOOKUP_TABLE default"<<"\n";
    for(int k = 0; k < num_of_elems; k++)
    {
        os<< block_index <<"\n";
    }

    os.close();
    return;


}


void MeshOpt::getFromInput(tbox::Pointer<tbox::Database> db){
}

void MeshOpt::putToDatabase(tbox::Pointer<tbox::Database> db)
{

}

#endif


