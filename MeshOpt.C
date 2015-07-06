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


/*************************************************************************
*    MeshOpt 类的构造函数.	 	                                *
*************************************************************************/
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

    d_flag = 0;


    // 为影像区的宽度赋值.
    d_zeroghosts = hier::IntVector<NDIM>(0);
    d_oneghosts  = hier::IntVector<NDIM>(1);

    // 创建所有变量及数据片索引号, 注册可视化数据片.
    registerModelVariables();

}

/*************************************************************************
*    MeshOpt 类的空析构函数.	 	                                *
*************************************************************************/

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
    d_coords = new pdat::NodeVariable<NDIM,double>("coordinates", NDIM);
    d_fixed = new pdat::NodeVariable<NDIM,bool>("fixed",1);

    // 当前值上下文, 新值上下文, 演算上下文.
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
                                                    d_zeroghosts);

    d_fixed_info_id =
            variable_db->registerVariableAndContext(d_fixed,
                                                    d_current,
                                                    d_zeroghosts);

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
        intc->registerInitPatchData(d_fixed_info_id);

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

    }else if(intc_name=="OPTIMIZE_MESH") { // 数值构件: 优化网格
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

/*************************************************************************
* 初值构件: 初始化网格结点坐标
*************************************************************************/
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

        tbox::Pointer< pdat::NodeData<NDIM,bool> > fixed_info
                = patch.getPatchData(d_fixed_info_id);

        // 赋初值
        fixed_info->fillAll(false);

        // 设定结点类型
        setNodeInfo(patch);
    }

}


/*************************************************************************
* 步长构件: 计算并返回网格片上的稳定时间步长.  
*************************************************************************/
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

void MeshOpt::computeOnPatch(hier::Patch<NDIM>& patch,
                             const double  time,
                             const double  dt,
                             const bool    initial_time,
                             const string& intc_name)
{

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(intc_name=="WRITE_MESH"|| "DISTRUB_MESH" || "OPTIMIZE_MESH");
#endif

    if(intc_name =="WRITE_MESH")
    {
        writeToVTK(patch,time,dt,initial_time);
    }
    else if(intc_name =="DISTRUB_MESH")
    {
        disturbMesh(patch,time,dt,initial_time);
    }
     else  if(intc_name =="OPTIMIZE_MESH")
   {
        trustRegion(patch);
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



/*
 ********************************************************************
 * 从数据库中读取信息或把信息存放到数据库中                                               *
 ********************************************************************
 */
void MeshOpt::getFromInput(tbox::Pointer<tbox::Database> db)
{

}

void MeshOpt::putToDatabase(tbox::Pointer<tbox::Database> db)
{

}


/********************************************************************
 * 网格处理                                             *
 ********************************************************************/
 void MeshOpt::setNodeInfo(hier::Patch<NDIM>& patch)
{
    const hier::Index<NDIM> ifirst=patch.getBox().lower();
    const hier::Index<NDIM> ilast =patch.getBox().upper();

    int cols = ilast(0) - ifirst(0)+2;
    int rows = ilast(1) - ifirst(1)+2;

     int num_of_nodes = cols*rows;

    tbox::Array<bool> d_fixed_info(num_of_nodes,true);

    for(int row = 0; row < rows; row++)
    {
        for(int col = 0; col < cols; col++)
        {
            if(row == 0 || row == rows-1 || col == 0 || col == cols-1)
            {
                d_fixed_info[row*cols+col] = true;
            }
            else
                d_fixed_info[row*cols+col] = false;
        }
    }

    tbox::Pointer< pdat::NodeData<NDIM,bool> > fixed_info
            = patch.getPatchData(d_fixed_info_id);

    int count=-1;
    for(pdat::NodeIterator<NDIM> ic((*fixed_info).getBox()); ic; ic++)
    {
       (*fixed_info)(ic(),0) = d_fixed_info[++count];
    }
}


MeshImpl * MeshOpt::createLocalMesh(hier::Patch<NDIM> & patch)
{
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


    tbox::Pointer< pdat::NodeData<NDIM,bool> > fixed_info
            = patch.getPatchData(d_fixed_info_id);

    MeshImpl * mesh = new  MeshImpl(num_of_nodes, num_of_elems,
                                    QUADRILATERAL,
                                    fixed_info.getPointer()->getPointer(),
                                    node.getPointer(),
                                    elem.getPointer());
    return mesh;
}

void MeshOpt::transformMeshtoPatch(MeshImpl * mesh, hier::Patch<NDIM>& patch, MsqError& err)
{
    std::vector<Mesh::VertexHandle>  vertices;
    mesh->get_all_vertices(vertices,err);
    size_t num_of_vertices = vertices.size();
//    std::cout << num_of_vertices << std::endl;

    std::vector<MsqVertex> coords(num_of_vertices);
    mesh->vertices_get_coordinates(arrptr(vertices),arrptr(coords),num_of_vertices,err);

    tbox::Pointer< pdat::NodeData<NDIM,double> > coords_current
            = patch.getPatchData(d_coords_current_id);

    int count = 0;
    for(pdat::NodeIterator<NDIM> ic((*coords_current).getBox()); ic; ic++)
    {
        (*coords_current)(ic(),0) = coords[count][0];
        (*coords_current)(ic(),1) =  coords[count][1];
        ++count;
    }

    return;

}


// 扰动网格
void MeshOpt::disturbMesh(hier::Patch<NDIM>& patch,
                          const double  time,
                          const double  dt,
                          const bool    initial_time)
{
    NULL_USE(dt);
    NULL_USE(time);
    NULL_USE(initial_time);

    tbox::Pointer< pdat::NodeData<NDIM,double> > coords_current
            = patch.getPatchData(d_coords_current_id);

    tbox::Pointer< pdat::NodeData<NDIM,bool> > fixed_info
            = patch.getPatchData(d_fixed_info_id);

    int count = -1;
    double dist = 0.01;
    for(pdat::NodeIterator<NDIM> ic((*coords_current).getBox()); ic; ic++)
    {
        dist *= -1;
        if((*fixed_info)(ic(),0) == false)
        {
            (*coords_current)(ic(),0) += dist;
            (*coords_current)(ic(),1) -= dist;
        }
        ++count;
    }

    // 表示已扰动
    d_flag = 1;

}

// 输出网格到 VTK 文件
void MeshOpt::writeToVTK(hier::Patch<NDIM>& patch,
                         const double  time,
                         const double  dt,
                         const bool    initial_time)
{
    NULL_USE(dt);
    NULL_USE(time);
    NULL_USE(initial_time);

    const tbox::Pointer< hier::BlockPatchGeometry<NDIM> > pgeom =
            patch.getPatchGeometry();

    int block_index = pgeom->getBlockNumber();

    int patch_index = patch.getPatchNumber();

    std::stringstream bi, pi, df;
    bi << block_index;
    pi << patch_index;
    df << d_flag;

    std::string file_name = df.str() + "_block_ " + bi.str()+ "_patch_" +  pi.str()  + ".vtk";


    MsqError err;
    MeshImpl * mesh = createLocalMesh(patch);
    mesh->write_vtk(file_name.c_str(), err);

    return;
}


void MeshOpt::trustRegion(hier::Patch<NDIM>& patch)
{
    MsqError err;
    MeshImpl* mesh = createLocalMesh(patch);

    PlanarDomain domain(PlanarDomain::XY);
    IdealWeightInverseMeanRatio inverse_mean_ratio(err);
    LPtoPTemplate obj_func(&inverse_mean_ratio,2,err);
    TrustRegion t_region(&obj_func);
    t_region.use_global_patch();
    TerminationCriterion tc_inner;
    tc_inner.add_absolute_gradient_L2_norm(1e-6);
    tc_inner.add_iteration_limit(1);
    t_region.set_inner_termination_criterion(&tc_inner);
    QualityAssessor m_ratio_qa(&inverse_mean_ratio);
    m_ratio_qa.disable_printing_results();

    InstructionQueue queue;
    queue.add_quality_assessor(&m_ratio_qa,err);
    queue.set_master_quality_improver(&t_region,err);
    queue.add_quality_assessor(&m_ratio_qa,err);
    MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(mesh,&domain);
    queue.run_instructions(&mesh_and_domain,err);
    tbox::pout<< "Shape optimization completed." << std::endl;

//    std::string file_name = "2__patch_3.vtk";
//    mesh->write_vtk(file_name.c_str(), err);

    transformMeshtoPatch(mesh,patch,err);



    d_flag =2;

    return;
}

#endif


