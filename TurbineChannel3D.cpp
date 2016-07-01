#include "TurbineChannel3D.h"
#include "lattice_vars.h"
#include <set>
#include <map>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <stdexcept>
#include <cmath>
#include <stdlib.h>


// for debugging
#include <iostream>

using namespace std;

const string TurbineChannel3D::params_file="params.lbm";
const string TurbineChannel3D::snl_file="snl.lbm";
const string TurbineChannel3D::inl_file="inl.lbm";
const string TurbineChannel3D::onl_file="onl.lbm";

TurbineChannel3D::TurbineChannel3D(const int rk, const int sz):
rank(rk), size(sz)
{

    const char * altDen = getenv("LBM_FN_DEN");
    const char * altSfx = getenv("LBM_FN_SFX");
    const char * altUxf = getenv("LBM_FN_UX");
    const char * altUyf = getenv("LBM_FN_UY");
    const char * altUzf = getenv("LBM_FN_UZ");
    // set defaults
    if (!altDen) { altDen = "density"; }
    if (!altSfx) { altSfx = ".b_dat"; }
    if (!altUxf) { altUxf = "ux"; }
    if (!altUyf) { altUyf = "uy"; }
    if (!altUzf) { altUzf = "uz"; }

    densityFileStub.assign(altDen);
    fileSuffix.assign(altSfx);
    ux_FileStub.assign(altUxf);
    uy_FileStub.assign(altUyf);
    uz_FileStub.assign(altUzf);

    tag_d = 666; tag_u = 999;
    read_input_file(params_file);
    initialize_lattice_data();
    initialize_local_partition_variables();
    initialize_mpi_buffers();
    vtk_ts = 0;
}

TurbineChannel3D::~TurbineChannel3D(){
    delete [] fEven;
    delete [] fOdd;
    delete [] inl;
    delete [] onl;
    delete [] snl;
    delete [] u_bc;
    
    delete [] ghost_in_m;
    delete [] ghost_out_m;
    delete [] ghost_in_p;
    delete [] ghost_out_p;
    
    delete [] rho_l;
    delete [] ux_l;
    delete [] uy_l;
    delete [] uz_l;
    
}

void TurbineChannel3D::write_bc_arrays(MPI_Comm comm){
  MPI_File fh_snl, fh_inl, fh_onl;
  MPI_Status mpi_s1, mpi_s2, mpi_s3;
  // open MPI data file
    MPI_File_open(comm,"snl.b_dat",
      MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&fh_snl);
    MPI_File_open(comm,"inl.b_dat",
      MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&fh_inl);
    MPI_File_open(comm,"onl.b_dat",
      MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&fh_onl);
  // write your chunk of the data
    int offset_s = firstSlice*Nx*Ny*sizeof(int);
    MPI_File_write_at(fh_snl,offset_s,snl+HALO*Nx*Ny,numEntries,MPI_INT,&mpi_s1);
    MPI_File_write_at(fh_inl,offset_s,inl+HALO*Nx*Ny,numEntries,MPI_INT,&mpi_s2);
    MPI_File_write_at(fh_onl,offset_s,onl+HALO*Nx*Ny,numEntries,MPI_INT,&mpi_s3);

  // close the MPI data file
    MPI_File_close(&fh_snl);
    MPI_File_close(&fh_inl);
    MPI_File_close(&fh_onl);
}

void TurbineChannel3D::write_data(MPI_Comm comm, bool isEven){
    
    stringstream ts_ind;
    string ts_ind_str;
    string density_fn, ux_fn, uy_fn, uz_fn;
    MPI_File fh_rho, fh_ux, fh_uy, fh_uz;
    MPI_Status mpi_s1, mpi_s2, mpi_s3,mpi_s4;
    
        
    const float * RESTRICT fIn;
    int streamNum;
    if (isEven){
        fIn = fEven;
        streamNum = 6;
    }else{
        fIn = fOdd;
        streamNum = 9;
    }
    
    int nnodes = this->nnodes;
    int numSpd = this->numSpd;
    int Nx = this->Nx;
    int Ny = this->Ny;
    const float * RESTRICT ex = this->ex;
    const float * RESTRICT ey = this->ey;
    const float * RESTRICT ez = this->ez;
    float * RESTRICT ux_l = this->ux_l;
    float * RESTRICT uy_l = this->uy_l;
    float * RESTRICT uz_l = this->uz_l;
    float * RESTRICT rho_l = this->rho_l;
    int totalSlices = this->totalSlices;
    const int * RESTRICT snl = this->snl;
    int numMySlices = this->numMySlices;
    int numEntries = Nx*Ny*numMySlices;
   
    
      for(int z = HALO;z<(totalSlices-HALO);z++){
          for(int y = 0;y<Ny;y++){
              for(int x = 0;x<Nx;x++){
                  int tid_l, tid_g;
                  float tmp_rho, tmp_ux, tmp_uy, tmp_uz;
                  tid_l = x+y*Nx+(z-HALO)*Nx*Ny;
                  tmp_rho = 0; tid_g = x+y*Nx+z*Nx*Ny;
                  tmp_ux = 0; tmp_uy = 0; tmp_uz = 0;
                  
                  for(int spd=0;spd<numSpd;spd++){
                      float f = fIn[getIdx(nnodes,numSpd,tid_g,spd)];
                      tmp_rho+=f;
                      tmp_ux+=ex[spd]*f;
                      tmp_uy+=ey[spd]*f;
                      tmp_uz+=ez[spd]*f;
                  }
                  rho_l[tid_l]=tmp_rho;
                  if(snl[tid_g]==1){
                      ux_l[tid_l]=0.; uy_l[tid_l]=0.; uz_l[tid_l]=0.;
                  }else{
                      ux_l[tid_l]=tmp_ux*(1./tmp_rho);
                      uy_l[tid_l]=tmp_uy*(1./tmp_rho);
                      uz_l[tid_l]=tmp_uz*(1./tmp_rho);
                  }
              }
          }
      }
      
      
      ts_ind << vtk_ts;
      density_fn = densityFileStub+ts_ind.str()+fileSuffix;
      ux_fn = ux_FileStub+ts_ind.str()+fileSuffix;
      uy_fn = uy_FileStub+ts_ind.str()+fileSuffix;
      uz_fn = uz_FileStub+ts_ind.str()+fileSuffix;
      
      // open MPI file for parallel IO
      MPI_File_open(comm,(char*)density_fn.c_str(),
      MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&fh_rho);
      
      MPI_File_open(comm,(char*)ux_fn.c_str(),
      MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&fh_ux);
      
      MPI_File_open(comm,(char*)uy_fn.c_str(),
      MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&fh_uy);
      
      MPI_File_open(comm,(char*)uz_fn.c_str(),
      MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&fh_uz);
      
      
      
      //write your chunk of data
      MPI_File_write_at(fh_rho,offset,rho_l,numEntries,MPI_FLOAT,&mpi_s1);
      MPI_File_write_at(fh_ux,offset,ux_l,numEntries,MPI_FLOAT,&mpi_s2);
      MPI_File_write_at(fh_uy,offset,uy_l,numEntries,MPI_FLOAT,&mpi_s3);
      MPI_File_write_at(fh_uz,offset,uz_l,numEntries,MPI_FLOAT,&mpi_s4);
      
      //close the files
      MPI_File_close(&fh_rho);
      MPI_File_close(&fh_ux);
      MPI_File_close(&fh_uy);
      MPI_File_close(&fh_uz);
      
      vtk_ts++; // increment the dump counter...
      
   }


    void TurbineChannel3D::D3Q15_process_slices(bool isEven, const int firstSlice, const int lastSlice, int streamNum, int waitNum){
    
    
    const float * RESTRICT fIn;
    float * RESTRICT fOut;
    int writeWaitNum;
    
    if(isEven){
        fIn = fEven; fOut = fOdd; writeWaitNum=6;
    }else{
        fIn = fOdd; fOut = fEven; writeWaitNum=9;
    }
    
    // local copies of class data members needed for acc compiler
    const int* RESTRICT inl = this->inl;
    const int* RESTRICT onl = this->onl;
    const int* RESTRICT snl = this->snl;
    const float* RESTRICT u_bc = this->u_bc;
    int Ny = this->Ny;
    int Nx = this->Nx;
    float omega_l = this->omega;
    int nnodes = this->nnodes;
    float rho_lbm = this->rho_lbm;
    

    float Cs = this->Cs;
    
   
   
    const int numSpd=15;
   
    for(int Z=firstSlice;Z<lastSlice;Z++){
        for(int Y=0;Y<Ny;Y++){
            for(int X=0;X<Nx;X++){
                float f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14;
                float cu,rho,ux,uy,uz,w;
                int X_t,Y_t,Z_t,tid_t,tid;
                
                tid=X+Y*Nx+Z*Nx*Ny;
                
                //load the data into registers
                f0=fIn[getIdx(nnodes, numSpd, tid,0)]; f1=fIn[getIdx(nnodes, numSpd, tid,1)];
                f2=fIn[getIdx(nnodes, numSpd, tid,2)]; f3=fIn[getIdx(nnodes, numSpd, tid,3)];
                f4=fIn[getIdx(nnodes, numSpd, tid,4)]; f5=fIn[getIdx(nnodes, numSpd, tid,5)];
                f6=fIn[getIdx(nnodes, numSpd, tid,6)]; f7=fIn[getIdx(nnodes, numSpd, tid,7)];
                f8=fIn[getIdx(nnodes, numSpd, tid,8)]; f9=fIn[getIdx(nnodes, numSpd, tid,9)];
                f10=fIn[getIdx(nnodes, numSpd, tid,10)]; f11=fIn[getIdx(nnodes, numSpd, tid,11)];
                f12=fIn[getIdx(nnodes, numSpd, tid,12)]; f13=fIn[getIdx(nnodes, numSpd, tid,13)];
                f14=fIn[getIdx(nnodes, numSpd, tid,14)];
                
                //compute density
                rho = f0+f1+f2+f3+f4+f5+f6+f7+f8+f9+f10+f11+f12+f13+f14;
                ux=f1-f2+f7-f8+f9-f10+f11-f12+f13-f14; ux/=rho;
                uy=f3-f4+f7+f8-f9-f10+f11+f12-f13-f14; uy/=rho;
                uz=f5-f6+f7+f8+f9+f10-f11-f12-f13-f14; uz/=rho;
                
                	// set macroscopic boundary conditions
		if(inl[tid]==1){
		  ux=0;uy=0; uz=u_bc[tid];
		  //set rho based on uz
		  rho = (1.F/(1.F-uz))*(2.0F*(f6+f11+f12+f13+f14)+(f0+f1+f2+f3+f4));

		}
		if(onl[tid]==1){
		  ux=0.; uy=0.; rho=rho_lbm;
		  uz = -1.F+((2.F*(f5+f7+f8+f9+f10)+(f0+f1+f2+f3+f4)))/rho;		 
		}
		if(snl[tid]==1){
		  ux=0.; uy=0.; uz=0.;
		}
	
	
		 //everyone compute equilibrium
		float fe0,fe1,fe2,fe3,fe4,fe5,fe6,fe7,fe8,fe9,fe10,fe11,fe12,fe13,fe14;
		//speed 0 ex=ey=ez=0 w=2./9.
	
		fe0=rho*(2.F/9.F)*(1.F-1.5F*(ux*ux+uy*uy+uz*uz));
	

		//speed 1 ex=1 ey=ez=0 w=1./9.
		cu=3.F*(1.F*ux);
		fe1=rho*(1.F/9.F)*(1.F+cu+0.5F*(cu*cu)-
				1.5F*(ux*ux+uy*uy+uz*uz));
	      

		//speed 2 ex=-1 ey=ez=0 w=1./9.
		cu=3.F*((-1.F)*ux);
		fe2=rho*(1.F/9.F)*(1.F+cu+0.5F*(cu*cu)-
				1.5F*(ux*ux+uy*uy+uz*uz));
	      

		//speed 3 ex=0 ey=1 ez=0 w=1./9.
		cu=3.F*(1.F*uy);
		fe3=rho*(1.F/9.F)*(1.F+cu+0.5F*(cu*cu)-
				1.5F*(ux*ux+uy*uy+uz*uz));
	      

		//speed 4 ex=0 ey=-1 ez=0 w=1./9.
		cu=3.F*(-1.F*uy);
		fe4=rho*(1.F/9.F)*(1.F+cu+0.5F*(cu*cu)-
				1.5F*(ux*ux+uy*uy+uz*uz));
	

		//speed 5 ex=ey=0 ez=1 w=1./9.
		cu=3.F*(1.F*uz);
		fe5=rho*(1.F/9.F)*(1.F+cu+0.5F*(cu*cu)-
				1.5F*(ux*ux+uy*uy+uz*uz));
	

		//speed 6 ex=ey=0 ez=-1 w=1./9.
		cu=3.F*(-1.F*uz);
		fe6=rho*(1.F/9.F)*(1.F+cu+0.5F*(cu*cu)-
				1.5F*(ux*ux+uy*uy+uz*uz));
	
		//speed 7 ex=ey=ez=1 w=1./72.
		cu=3.F*(ux+uy+uz);
		fe7=rho*(1.F/72.F)*(1.F+cu+0.5F*(cu*cu)-
				  1.5F*(ux*ux+uy*uy+uz*uz));
	      

		//speed 8 ex=-1 ey=ez=1 w=1./72.
		cu=3.F*(-ux+uy+uz);
		fe8=rho*(1.F/72.F)*(1.F+cu+0.5F*(cu*cu)-
				  1.5F*(ux*ux+uy*uy+uz*uz));
	

		//speed 9 ex=1 ey=-1 ez=1 w=1./72.
		cu=3.F*(ux-uy+uz);
		fe9=rho*(1.F/72.F)*(1.F+cu+0.5F*(cu*cu)-
				  1.5F*(ux*ux+uy*uy+uz*uz));
	

		//speed 10 ex=-1 ey=-1 ez=1 w=1/72
		cu=3.F*(-ux-uy+uz);
		fe10=rho*(1.F/72.F)*(1.F+cu+0.5F*(cu*cu)-
				  1.5F*(ux*ux+uy*uy+uz*uz));
	      

		//speed 11 ex=1 ey=1 ez=-1 w=1/72
		cu=3.F*(ux+uy-uz);
		fe11=rho*(1.F/72.F)*(1.F+cu+0.5F*(cu*cu)-
				  1.5F*(ux*ux+uy*uy+uz*uz));
	

		//speed 12 ex=-1 ey=1 ez=-1 w=1/72
		cu=3.F*(-ux+uy-uz);
		fe12=rho*(1.F/72.F)*(1.F+cu+0.5F*(cu*cu)-
				  1.5F*(ux*ux+uy*uy+uz*uz));
	      

		//speed 13 ex=1 ey=ez=-1 w=1/72
		cu=3.F*(ux-uy-uz);
		fe13=rho*(1.F/72.F)*(1.F+cu+0.5F*(cu*cu)-
				  1.5F*(ux*ux+uy*uy+uz*uz));
	      

		//speed 14 ex=ey=ez=-1 w=1/72
		cu=3.F*(-ux-uy-uz);
		fe14=rho*(1.F/72.F)*(1.F+cu+0.5F*(cu*cu)-
				  1.5F*(ux*ux+uy*uy+uz*uz));
	
		// if on inlet or outlet, compute and bounce-back non-equilibrium part of f.
		// see referenecs on regularized boundary conditions for full details on theory.
		if((inl[tid]==1)|(onl[tid]==1)){

		  float ft1,ft2,ft3,ft4,ft5,ft6,ft7,ft8,ft9,ft10,ft11,ft12,ft13,ft14;
		  if(inl[tid]==1){
		    //adjust fIn for the unknown velocities: 5,7,8,9,10
		    //bounce-back of non-equilibrium parts
		    //f5, bb_spd=f6
		    f5=fe5+(f6-fe6); //fIn[5*nnodes+tid]=f5;
		    //f7, bb_spd=f14
		    f7=fe7+(f14-fe14); //fIn[7*nnodes+tid]=f7;
		    //f8, bb_spd=f13
		    f8=fe8+(f13-fe13); //fIn[8*nnodes+tid]=f8;
		    //f9, bb_spd=f12
		    f9=fe9+(f12-fe12); //fIn[9*nnodes+tid]=f9;
		    //f10, bb_spd=f11
		    f10=fe10+(f11-fe11); //fIn[10*nnodes+tid]=f10;
		  }else{
		     f6=fe6+(f5-fe5); 
		    f11=fe11+(f10-fe10); 
		    f12=fe12+(f9-fe9); 
		    f13=fe13+(f8-fe8); 
		    f14=fe14+(f7-fe7); 
		  }
		  //ft0=f0-fe0;
		  ft1=f1-fe1; 
		  ft2=f2-fe2;
		  ft3=f3-fe3;
		  ft4=f4-fe4;
		  ft5=f5-fe5;
		  ft6=f6-fe6;
		  ft7=f7-fe7;
		  ft8=f8-fe8;
		  ft9=f9-fe9;
		  ft10=f10-fe10;
		  ft11=f11-fe11;
		  ft12=f12-fe12;
		  ft13=f13-fe13;
		  ft14=f14-fe14;

		  //now, multiply by f# = ((ft#)*Q_flat)*Q_flat'
		  f0= - ft1/3.F - ft2/3.F - ft3/3.F - ft4/3.F - ft5/3.F - ft6/3.F - ft7 - ft8 - ft9 - ft10 - ft11 - ft12 - ft13 - ft14; 
		  f1=(2.F*ft1)/3.F + (2.F*ft2)/3.F - ft3/3.F - ft4/3.F - ft5/3.F - ft6/3.F; 
		  f2=(2.F*ft1)/3.F + (2.F*ft2)/3.F - ft3/3.F - ft4/3.F - ft5/3.F - ft6/3.F; 
		  f3=(2.F*ft3)/3.F - ft2/3.F - ft1/3.F + (2.F*ft4)/3.F - ft5/3.F - ft6/3.F; 
		  f4=(2.F*ft3)/3.F - ft2/3.F - ft1/3.F + (2.F*ft4)/3.F - ft5/3.F - ft6/3.F; 
		  f5=(2.F*ft5)/3.F - ft2/3.F - ft3/3.F - ft4/3.F - ft1/3.F + (2.F*ft6)/3.F; 
		  f6=(2.F*ft5)/3.F - ft2/3.F - ft3/3.F - ft4/3.F - ft1/3.F + (2.F*ft6)/3.F; 
		  f7=(2.F*ft1)/3.F + (2.F*ft2)/3.F + (2.F*ft3)/3.F + (2.F*ft4)/3.F + (2.F*ft5)/3.F + (2.F*ft6)/3.F + 8.F*ft7 + 8.F*ft14;
		  f8= (2.F*ft1)/3.F + (2.F*ft2)/3.F + (2.F*ft3)/3.F + (2.F*ft4)/3.F + (2.F*ft5)/3.F + (2.F*ft6)/3.F + 8.F*ft8 + 8.F*ft13;
		  f9= (2.F*ft1)/3.F + (2.F*ft2)/3.F + (2.F*ft3)/3.F + (2.F*ft4)/3.F + (2.F*ft5)/3.F + (2.F*ft6)/3.F + 8.F*ft9 + 8.F*ft12;
		  f10= (2.F*ft1)/3.F + (2.F*ft2)/3.F + (2.F*ft3)/3.F + (2.F*ft4)/3.F + (2.F*ft5)/3.F + (2.F*ft6)/3.F + 8.F*ft10 + 8.F*ft11;
		  f11= (2.F*ft1)/3.F + (2.F*ft2)/3.F + (2.F*ft3)/3.F + (2.F*ft4)/3.F + (2.F*ft5)/3.F + (2.F*ft6)/3.F + 8.F*ft10 + 8.F*ft11;
		  f12= (2.F*ft1)/3.F + (2.F*ft2)/3.F + (2.F*ft3)/3.F + (2.F*ft4)/3.F + (2.F*ft5)/3.F + (2.F*ft6)/3.F + 8.F*ft9 + 8.F*ft12;
		  f13= (2.F*ft1)/3.F + (2.F*ft2)/3.F + (2.F*ft3)/3.F + (2.F*ft4)/3.F + (2.F*ft5)/3.F + (2.F*ft6)/3.F + 8.F*ft8 + 8.F*ft13;
		  f14= (2.F*ft1)/3.F + (2.F*ft2)/3.F + (2.F*ft3)/3.F + (2.F*ft4)/3.F + (2.F*ft5)/3.F + (2.F*ft6)/3.F + 8.F*ft7 + 8.F*ft14;

		  //update fIn for all velocities based on strain tensor
		  //f0, still equals 0..
		  cu = 9.F/2.F; w = 1.F/9.F;

		  //fIn[..] = fe#+f#
		  f0=fe0+f0;

		  f1=fe1+f1*(cu)*w;
		  f2=fe2+f2*(cu)*w;
		  f3=fe3+f3*cu*w;
		  f4=fe4+f4*cu*w;
		  f5=fe5+f5*cu*w;
		  f6=fe6+f6*cu*w;
		  w = 1.F/72.F;
		  f7=fe7+f7*cu*w;
		  f8=fe8+f8*cu*w;
		  f9=fe9+f9*cu*w;
		  f10=fe10+f10*cu*w;
		  f11=fe11+f11*cu*w;
		  f12=fe12+f12*cu*w;
		  f13=fe13+f13*cu*w;
		  f14=fe14+f14*cu*w;
		}

                // compute strain tensor and update local viscosity/relaxation factor accordingly
                float s11,s12,s13,s22,s23,s33;
		s11=(f1-fe1)+(f2-fe2)+(f7-fe7)+(f8-fe8)+(f9-fe9)+(f10-fe10)+(f11-fe11)+(f12-fe12)+(f13-fe13)+(f14-fe14);
		s12=(f7-fe7)-(f8-fe8)-(f9-fe9)+(f10-fe10)+(f11-fe11)-(f12-fe12)-(f13-fe13)+(f14-fe14);
		s13=(f7-fe7)-(f8-fe8)+(f9-fe9)-(f10-fe10)-(f11-fe11)+(f12-fe12)-(f13-fe13)+(f14-fe14);
		s22=(f3-fe3)+(f4-fe4)+(f7-fe7)+(f8-fe8)+(f9-fe9)+(f10-fe10)+(f11-fe11)+(f12-fe12)+(f13-fe13)+(f14-fe14);
		s23=(f7-fe7)+(f8-fe8)-(f9-fe9)-(f10-fe10)-(f11-fe11)-(f12-fe12)+(f13-fe13)+(f14-fe14);
		s33=(f5-fe5)+(f6-fe6)+(f7-fe7)+(f8-fe8)+(f9-fe9)+(f10-fe10)+(f11-fe11)+(f12-fe12)+(f13-fe13)+(f14-fe14);
		
                   
		float nu = 1.F/omega_l; nu-=0.5F; nu/=3.F;
		float P = s11*s11+2.F*s12*s12+2.F*s13*s13+s22*s22+2.F*s23*s23+s33*s33;
		P = sqrt(P);
		P*=Cs;
		P+=nu*nu;
		P=sqrt(P);
		P-=nu;

		//compute turbulent nu
		float nu_e = (1.F/(6.F))*P;
		//update omega
		    
		float omega = 1.F/(3.F*(nu+nu_e)+0.5F); //<-- shadows class data member this->omega

		//everyone (except solid nodes) relax...
                if(snl[tid]!= 1){
		  f0=f0-omega*(f0-fe0);
		  f1=f1-omega*(f1-fe1);
		  f2=f2-omega*(f2-fe2);
		  f3=f3-omega*(f3-fe3);
		  f4=f4-omega*(f4-fe4);
		  f5=f5-omega*(f5-fe5);
		  f6=f6-omega*(f6-fe6);
		  f7=f7-omega*(f7-fe7);
		  f8=f8-omega*(f8-fe8);
		  f9=f9-omega*(f9-fe9);
		  f10=f10-omega*(f10-fe10);
		  f11=f11-omega*(f11-fe11);
		  f12=f12-omega*(f12-fe12);
		  f13=f13-omega*(f13-fe13);
		  f14=f14-omega*(f14-fe14);
                } else {
                    // snl bounce-back
                    cu = f1; f1 = f2; f2 = cu;
                    cu = f3; f3 = f4; f4 = cu;
                    cu = f5; f5 = f6; f6 = cu;
                    cu = f7; f7 = f14; f14 = cu;
                    cu = f8; f8 = f13; f13 = cu;
                    cu = f9; f9 = f12; f12 = cu;
                    cu = f10; f10 = f11; f11 = cu; 
                }

		// stream data
                
                //speed 0 ex=ey=ez=0
                //fOut[tid]=f0;
                fOut[getIdx(nnodes, numSpd, tid,0)]=f0;
                
                //speed 1 ex=1 ey=ez=0
                X_t=X+1; Y_t=Y; Z_t=Z;
                if(X_t==Nx) X_t=0;
                tid_t=X_t+Y_t*Nx+Z_t*Nx*Ny;
                //	fOut[Nx*Ny*Nz+tid_t]=f1;
                fOut[getIdx(nnodes, numSpd, tid_t,1)]=f1;
                
                //speed 2 ex=-1 ey=ez=0;
                X_t=X-1; Y_t=Y; Z_t=Z;
                if(X_t<0) X_t=(Nx-1);
                tid_t=X_t+Y_t*Nx+Z_t*Nx*Ny;
                //	fOut[2*Nx*Ny*Nz+tid_t]=f2;
                fOut[getIdx(nnodes, numSpd, tid_t,2)]=f2;
                
                //speed 3 ex=0 ey=1 ez=0
                X_t=X; Y_t=Y+1; Z_t=Z;
                if(Y_t==Ny) Y_t=0;
                tid_t=X_t+Y_t*Nx+Z_t*Nx*Ny;
                //	fOut[3*Nx*Ny*Nz+tid_t]=f3;
                fOut[getIdx(nnodes, numSpd, tid_t,3)]=f3;
                
                //speed 4 ex=0 ey=-1 ez=0
                X_t=X; Y_t=Y-1; Z_t=Z;
                if(Y_t<0) Y_t=(Ny-1);
                tid_t=X_t+Y_t*Nx+Z_t*Nx*Ny;
                ///	fOut[4*Nx*Ny*Nz+tid_t]=f4;
                fOut[getIdx(nnodes, numSpd, tid_t,4)]=f4;
                
                
                //speed 5 ex=ey=0 ez=1
                X_t=X; Y_t=Y; Z_t=Z+1;
                //	if(Z_t==Nz) Z_t=0;
                tid_t=X_t+Y_t*Nx+Z_t*Nx*Ny;
                //fOut[5*Nx*Ny*Nz+tid_t]=f5;
                fOut[getIdx(nnodes, numSpd, tid_t,5)]=f5;
                
                //speed 6 ex=ey=0 ez=-1
                X_t=X; Y_t=Y; Z_t=Z-1;
                //	if(Z_t<0) Z_t=(Nz-1);
                tid_t=X_t+Y_t*Nx+Z_t*Nx*Ny;
                //	fOut[6*Nx*Ny*Nz+tid_t]=f6;
                fOut[getIdx(nnodes, numSpd, tid_t,6)]=f6;
                
                //speed 7 ex=ey=ez=1
                X_t=X+1; Y_t=Y+1; Z_t=Z+1;
                if(X_t==Nx) X_t=0;
                if(Y_t==Ny) Y_t=0;
                //	if(Z_t==Nz) Z_t=0;
                tid_t=X_t+Y_t*Nx+Z_t*Nx*Ny;
                //	fOut[7*Nx*Ny*Nz+tid_t]=f7;
                fOut[getIdx(nnodes, numSpd, tid_t,7)]=f7;
                
                //speed 8 ex=-1 ey=1 ez=1
                X_t=X-1; Y_t=Y+1; Z_t=Z+1;
                if(X_t<0) X_t=(Nx-1);
                if(Y_t==Ny) Y_t=0;
                //	if(Z_t==Nz) Z_t=0;
                tid_t=X_t+Y_t*Nx+Z_t*Nx*Ny;
                //	fOut[8*Nx*Ny*Nz+tid_t]=f8;
                fOut[getIdx(nnodes, numSpd, tid_t,8)]=f8;
                
                //speed 9 ex=1 ey=-1 ez=1
                X_t=X+1; Y_t=Y-1; Z_t=Z+1;
                if(X_t==Nx) X_t=0;
                if(Y_t<0) Y_t=(Ny-1);
                //	if(Z_t==Nz) Z_t=0;
                tid_t=X_t+Y_t*Nx+Z_t*Nx*Ny;
                //	fOut[9*Nx*Ny*Nz+tid_t]=f9;
                fOut[getIdx(nnodes, numSpd, tid_t,9)]=f9;
                
                //speed 10 ex=-1 ey=-1 ez=1
                X_t=X-1; Y_t=Y-1; Z_t=Z+1;
                if(X_t<0) X_t=(Nx-1);
                if(Y_t<0) Y_t=(Ny-1);
                //	if(Z_t==Nz) Z_t=0;
                tid_t=X_t+Y_t*Nx+Z_t*Nx*Ny;
                //	fOut[10*Nx*Ny*Nz+tid_t]=f10;
                fOut[getIdx(nnodes, numSpd, tid_t,10)]=f10;
                
                //speed 11 ex=1 ey=1 ez=-1
                X_t=X+1; Y_t=Y+1; Z_t=Z-1;
                if(X_t==Nx) X_t=0;
                if(Y_t==Ny) Y_t=0;
                //	if(Z_t<0) Z_t=(Nz-1);
                tid_t=X_t+Y_t*Nx+Z_t*Nx*Ny;
                //	fOut[11*Nx*Ny*Nz+tid_t]=f11;
                fOut[getIdx(nnodes, numSpd, tid_t,11)]=f11;
                
                //speed 12 ex=-1 ey=1 ez=-1
                X_t=X-1; Y_t=Y+1; Z_t=Z-1;
                if(X_t<0) X_t=(Nx-1);
                if(Y_t==Ny) Y_t=0;
                //	if(Z_t<0) Z_t=(Nz-1);
                tid_t=X_t+Y_t*Nx+Z_t*Nx*Ny;
                //	fOut[12*Nx*Ny*Nz+tid_t]=f12;
                fOut[getIdx(nnodes, numSpd, tid_t,12)]=f12;
                
                //speed 13 ex=1 ey=-1 ez=-1
                X_t=X+1; Y_t=Y-1; Z_t=Z-1;
                if(X_t==Nx) X_t=0;
                if(Y_t<0) Y_t=(Ny-1);
                //	if(Z_t<0) Z_t=(Nz-1);
                tid_t=X_t+Y_t*Nx+Z_t*Nx*Ny;
                //	fOut[13*Nx*Ny*Nz+tid_t]=f13;
                fOut[getIdx(nnodes, numSpd, tid_t,13)]=f13;
                
                //speed 14 ex=ey=ez=-1
                X_t=X-1; Y_t=Y-1; Z_t=Z-1;
                if(X_t<0) X_t=(Nx-1);
                if(Y_t<0) Y_t=(Ny-1);
                //	if(Z_t<0) Z_t=(Nz-1);
                tid_t=X_t+Y_t*Nx+Z_t*Nx*Ny;
                //	fOut[14*Nx*Ny*Nz+tid_t]=f14;
                
                //fOut[tid_t*numSpd+14]=f14;
                fOut[getIdx(nnodes, numSpd, tid_t,14)]=f14;
            }
        }
    }
}

    void TurbineChannel3D::stream_out_collect(bool isEven,const int z_start,float * RESTRICT buff_out, const int numStreamSpeeds, const int * RESTRICT  streamSpeeds, int streamNum){
    int Ny = this->Ny;
    int Nx = this->Nx;
    int numSpd = this->numSpd;
    int nnodes = this->nnodes;
    
    
    float * RESTRICT fIn_b;
    if(isEven) {
      fIn_b = fOdd;
    }else{
      fIn_b = fEven;
    }
    
    
    for(int z=0;z<HALO;z++){
        for(int y=0;y<Ny;y++){
            for(int x=0;x<Nx;x++){
                for(int spd=0;spd<numStreamSpeeds;spd++){
                    int tid_l = x+y*Nx+z*Nx*Ny; int tid_g = x+y*Nx+(z+z_start)*Nx*Ny;
                    int stream_spd=streamSpeeds[spd];
                    buff_out[tid_l*numStreamSpeeds+spd]=fIn_b[getIdx(nnodes, numSpd, tid_g,stream_spd)];
                }
            }
        }
    }


}

    void TurbineChannel3D::stream_in_distribute(bool isEven,const int z_start, const float * RESTRICT buff_in, const int numStreamSpeeds, const int * RESTRICT streamSpeeds, int streamNum){
    int Nx = this->Nx;
    int Ny = this->Ny;
    int numSpd = this->numSpd;
    int nnodes = this->nnodes;
    
    float * RESTRICT fOut_b;
    
    if (isEven) {
      fOut_b = fOdd;
    }else{
      fOut_b = fEven;
    }
    
   
    for(int z=0;z<HALO;z++){
        for(int y=0;y<Ny;y++){
            for(int x=0;x<Nx;x++){
                for(int spd=0;spd<numStreamSpeeds;spd++){
                    int tid_l=x+y*Nx+z*Nx*Ny; int tid_g = x+y*Nx+(z+z_start)*Nx*Ny;
                    int stream_spd=streamSpeeds[spd];
                    fOut_b[getIdx(nnodes, numSpd, tid_g,stream_spd)]=buff_in[tid_l*numStreamSpeeds+spd];
                }
            }
        }
    }
}

void TurbineChannel3D::take_lbm_timestep(bool isEven, MPI_Comm comm){
     // --- Start Sequential Dependency A -------
    int waitNum;
    if (isEven) {
       waitNum = 3;
    }else{
       waitNum = 2;
    }
    // collide and stream lower boundary slices
    D3Q15_process_slices(isEven,HALO,HALO+1,0,waitNum);
    
    // collect data from lower HALO slice z = 0
    stream_out_collect(isEven,0,ghost_out_m,numMspeeds,Mspeeds,0);
    
   
    // ---  Sequential Dependency A ------
    
    // --- Start Sequential Dependency B -----
    // collide and stream upper boundary slices
    D3Q15_process_slices(isEven,totalSlices-2*HALO,totalSlices-HALO,1,waitNum);
    
    // collect data from upper HALO slice; z = totalSlices-1
    stream_out_collect(isEven,totalSlices-HALO,ghost_out_p,numPspeeds,Pspeeds,1);
    
    // --- Start Sequential Dependency C -----
    int streamNum;
    if (isEven) {
       streamNum = 2;
    }else{
       streamNum = 3;
    }

    // collide and stream interior lattice points
    D3Q15_process_slices(isEven,HALO+1,totalSlices-2*HALO,streamNum,waitNum);


 // begin communication to ghost_p_in
    MPI_Isend(ghost_out_m,numHALO,MPI_FLOAT,nd_m,tag_d,comm, &rq_out1);
    MPI_Irecv(ghost_in_p,numHALO,MPI_FLOAT,nd_p,tag_d,comm,&rq_in1);
    // begin communication to ghost_m_in
    MPI_Isend(ghost_out_p,numHALO,MPI_FLOAT,nd_p,tag_u,comm,&rq_out2);
    MPI_Irecv(ghost_in_m,numHALO,MPI_FLOAT,nd_m,tag_u,comm,&rq_in2);
    // ---  Sequential Dependency B -----
    
    // --- End Sequential Dependency C -----
    // ensure communication of boundary lattice points is complete
   
    MPI_Wait(&rq_in1,&stat); // Sequential Dependency A
    MPI_Wait(&rq_in2,&stat); // Sequential Dependency B
   
    
    // copy data from i+1 partition into upper boundary slice
    // Sequential Dependency A  
    stream_in_distribute(isEven,totalSlices-2*HALO,ghost_in_p,numMspeeds,Mspeeds,0);
    
    // Sequential Dependency B
    // copy data from i-1 partition into lower boundary slice
    stream_in_distribute(isEven,HALO,ghost_in_m,numPspeeds,Pspeeds,1);
   
}

void TurbineChannel3D::initialize_mpi_buffers(){
    // integer for arithmetic into pointer arrays for nodes I need to
    // communicate to neighbors.
  firstNdm = Nx*Ny*HALO; // first node index on the lower boundary
  lastNdm = Nx*Ny*(HALO+1); // last node index on the lower boundary
  firstNdp = nnodes-2*(Nx*Ny*HALO); //first node index on the upper boundary 
  lastNdp = nnodes-(Nx*Ny*HALO); // last node index on the upper boundary
    
    
    
  numHALO = (Nx*Ny*numPspeeds*HALO); // number of nodes in each halo region
    
  ghost_in_m = new float[Nx*Ny*numMspeeds*HALO]; // buffer to hold incoming halo data from m neighbor
  ghost_out_m = new float[Nx*Ny*numMspeeds*HALO]; // buffer to hold outgoing halo data destined to m neighbor.
  ghost_in_p = new float[Nx*Ny*numPspeeds*HALO]; // buffer to hold incoming halo data from p neighbor
  ghost_out_p = new float[Nx*Ny*numPspeeds*HALO]; // buffer to hold outgoing halo data destined to p neighbor.
  offset = firstSlice*Nx*Ny*sizeof(float); // offset for writing into data file
  numEntries = numMySlices*Nx*Ny; // size of entry for data file
    
  rho_l = new float[numEntries]; // buffer to hold output data destined for data file
  ux_l = new float[numEntries]; // ditto
  uy_l = new float[numEntries]; // ditto
  uz_l = new float[numEntries]; // ditto
}

void TurbineChannel3D::initialize_local_partition_variables(){
    int tstRank = 0; // rank where I am currently testing
  // compute number of slices for this patition
    numMySlices = Nz/size;
    if(rank<(Nz%size))
    numMySlices+=1;
    
// identify first and last slice and count the total number of nodes (including HALO)
    firstSlice=(Nz/size)*rank;
    if((Nz%size)<rank){
        firstSlice+=(Nz%size);
    }else{
        firstSlice+=rank;
    }
    lastSlice = firstSlice+numMySlices-1;
    totalSlices=numMySlices+2*HALO;// add 2 HALO slices (HALO=1)
    nnodes = totalSlices*Nx*Ny;


    

    // allocate memory for dependent variables and BC arrays
    fEven = new float[nnodes*numSpd];
    fOdd = new float[nnodes*numSpd];
    snl = new int[nnodes];
    inl = new int[nnodes];
    onl = new int[nnodes];
    u_bc = new float[nnodes];

    // identify node to left (nd_m) and node to right (nd_p)
    nd_m = rank-1;
    if(nd_m<0)
    nd_m=(size-1);
    
    nd_p = rank+1;
    if(nd_p==size)
    nd_p=0;
    


    // initialize dependent variable arrays
    int tid;
    for(int z=0; z<totalSlices;z++){
        for(int y=0;y<Ny;y++){
            for(int x=0;x<Nx;x++){
                tid=x+y*Nx+z*Nx*Ny;
                for(int spd=0;spd<numSpd;spd++){
                    fEven[getIdx(nnodes, numSpd, tid,spd)]=rho_lbm*w[spd];
                    fOdd[getIdx(nnodes,numSpd,tid,spd)]=rho_lbm*w[spd];
                }
            }
        }
    }
    
    // populate myGlobalNodes set and globalToLocal map
    int l_id; //local id    
    for(int z = firstSlice;z<=lastSlice;z++){
      for(int y = 0;y<Ny;y++){
        for(int x = 0;x<Nx;x++){
          tid = x+y*Nx+z*Nx*Ny; // Global node number
          l_id = ((z-firstSlice)*Nx*Ny+x+y*Nx)+Nx*Ny*HALO; // node number accounting for HALO

          
          myGlobalNodes.insert(tid);
          globalToLocal[tid]=l_id;
           
        }
      }
    }


    // initialize all boundary condition node lists to zero...
    for(int nd=0;nd<nnodes;nd++){
        inl[nd]=0;
        onl[nd]=0;
        snl[nd]=0;
    }
    
    //inl, onl, snl and u_bc all need to be set based on data given in 
    //the respecive *.lbm file.  
   
   
    int numBCnode;
    int bcNode;
   

    
    // set snl
    fstream bcFile(snl_file.c_str(),ios::in);
    bcFile >> numBCnode;
    for(int i=0;i<numBCnode;i++){
      bcFile >> bcNode;

      // determine if the solid node is on this partition.
      if(myGlobalNodes.find(bcNode) != myGlobalNodes.end()){
         snl[globalToLocal[bcNode]] = 1;
      }
    }
    bcFile.close();

  

    // set inl 
    bcFile.open(inl_file.c_str(),ios::in);
    bcFile >> numBCnode;
    for(int i=0;i<numBCnode;i++){
      bcFile >> bcNode;

      // determine if the inlet node is on this partition.
      if(myGlobalNodes.find(bcNode) != myGlobalNodes.end()){
         inl[globalToLocal[bcNode]] = 1;
      }
    }
    bcFile.close();

    // set onl
    bcFile.open(onl_file.c_str(),ios::in);
    bcFile >> numBCnode;
    for(int i=0;i<numBCnode;i++){
      bcFile >> bcNode;

      // determine if the outlet node is on this partition.
      if(myGlobalNodes.find(bcNode) != myGlobalNodes.end()){
         onl[globalToLocal[bcNode]] = 1;
      }
    }
    bcFile.close();

    
    // initialize u_bc
   
    for(int z=0;z<totalSlices;z++){
        for(int y=0;y<Ny;y++){
            for(int x=0;x<Nx;x++){
                tid=x+y*Nx+z*Nx*Ny;
                if((inl[tid]==1)){
                     u_bc[tid]=umax_lbm; // constant inlet velocity
                }else{
                    u_bc[tid]=0.;
                }
            }
        }
    }
    
    // initialize outlet density
    rho_out = rho_lbm;
}

void TurbineChannel3D::read_input_file(const string input_file){
    ifstream input_params(input_file.c_str(),ios::in);
    if(!input_params.is_open())
    throw std::runtime_error("Could not open params file!");
    input_params >> LatticeType;
    input_params >> Num_ts;
    input_params >> ts_rep_freq;
    input_params >> Warmup_ts;
    input_params >> plot_freq;
    input_params >> Cs;
    input_params >> rho_lbm;
    input_params >> umax_lbm;
    input_params >> omega;
    input_params >> Nx;
    input_params >> Ny;
    input_params >> Nz;
    input_params.close();
}

void TurbineChannel3D::initialize_lattice_data(){
    switch(LatticeType){
        case(1):
        
        ex = ex15;
        ey = ey15;
        ez = ez15;
        w = w15;
        numPspeeds = numPspeedsD3Q15;
        numMspeeds = numMspeedsD3Q15;
        Mspeeds = MspeedsD3Q15;
        Pspeeds = PspeedsD3Q15; 
    }
}
