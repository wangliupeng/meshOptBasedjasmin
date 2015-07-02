##
## 文件名:	Makefile.in
## 软件包:	JASMIN applications
## 版权  :	(c) 2004-2010 北京应用物理与计算数学研究所
## 版本号:	$Revision: 170 $
## 修改  :	$Date: 2007-06-27 08:44:22  $
## 描述  :	makefile for Euler gas dynamics sample application
##

JASMIN	      = /usr/local/jasmin/3.0-dev

FORTRAN       = fortran
M4DIRS        = -DFORTDIR=$(FORTRAN)/$(NDIM)d -DJASMIN_FORTDIR=$(JASMIN)/include

          

default:        main-2d

include $(JASMIN)/config/Makefile.config

CPPFLAGS_EXTRA= -DNDIM=$(NDIM)

CXX_OBJS      = main_MeshOpt.o  MeshOpt.o MeshOptLevelIntegrator.o


main-2d:
		if test -f stamp-3d; then $(MAKE) clean; fi
		touch stamp-2d
		$(MAKE) NDIM=2 main2d

main2d:		$(CXX_OBJS) $(F2D_OBJS) $(LIBJASMINDEPEND)
		$(CXX) $(CXXFLAGS) $(CXXLD_FLAGS) $(CXX_OBJS) $(F2D_OBJS)	\
		$(LIBJASMIN2D) $(LIBJASMIN) $(CXXLDLIBS) -o main2d

clean:
		$(RM) *.f *.o main2d

redo:
		$(RM) main2d core *.ii *.int.c
		$(RM) -r ti_files ii_files

include Makefile.depend
