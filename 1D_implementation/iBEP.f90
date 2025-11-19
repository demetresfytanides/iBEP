           program column
!     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     +    modified by                                                 +
!     +    Dimitrios K. Fytanidis, Argonne National Laboratory         +
!     +                            email:dfytanidis@anl.gov            + 
!     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     +    by A. Martilli,     CIEMAT  SP 2040 MADRID                  +
!     +                        phone: ++34-9-13-46-62-99               +
!     +                        email:alberto.martilli@ciemat.es        +
!     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
           implicit none
           integer nzm
           parameter (nzm=128)
! flow variables
           real vx(nzm) ! x component of the wind
           real vy(nzm) ! y component of the wind
           real tke(nzm) ! tke
           real trc(nzm) ! trc
           real dpdx,dpdy ! external pressure gradient in x and y direction
           real utau_x,utau_y ! utau 
! variable needed in the calculation of the TKE and vertical diffusion  
           real dlk(nzm) ! mixing length
           real dls(nzm) ! dissipation length scale
           real cdz(nzm+1) ! vertical diffusion coefficients for momentum
	   real cdt(nzm+1) ! vertical diffusion coefficients for scalars
           real sh(nzm) ! shear production term
           real ceps,ck,cmu ! constants for the k-l scheme
           real pr ! prandtl number
           parameter (ceps=1/1.4,ck=0.4,cmu=0.09)
! variables needed in the resolution of the diffusion equation: here source/sinks for a generic variable c are defined as c*srim_c+srex_c
           real srim_vx(nzm) ! implicit part of the source/sinks terms of v1
           real srim_vy(nzm) ! implicit part of the source/sinks terms of v2
           real srim_tke(nzm) ! implicit part of the source/sinks terms of tke
           real srim_trc(nzm) ! implicit part of the source/sinks terms of tracer
           real srex_vx(nzm) ! implicit part of the source/sinks terms of v1
           real srex_vy(nzm) ! implicit part of the source/sinks terms of v2
           real srex_tke(nzm) ! implicit part of the source/sinks terms of tke
           real srex_trc(nzm) ! implicit part of the source/sinks terms of tracer
           real vl(nzm) ! fraction of air in each cell
           real sf(nzm+1) ! fraction of air at the interface between cells
! dispersive flux
           real dwtrc(nzm+1)
	 real duw(nzm+1)
           real aaa(4,4),bbb(4)

! drag coefficient
           real cdragx(nzm),cdragy(nzm) ! drag coefficient in x and y direction           
! building information
           real ss(nzm+1) ! ss(iz)=probability to have a building of height z(iz)
           real pb(nzm+1) ! pb(iz)=probbaility to have a building taller or equal to z(iz)
           real wx,wy ! distance between buildings at street level in the x and y direction respectively
           real bx,by ! building dimension in the x and y direction respectively
           real hmean ! mean building height
! emission at surface
           real emi
! numerical parameters
           real dz ! vertical resolution
           integer nz ! number of vertical levels
           real z(nzm+1) ! height of the faces of the numerical grid
           real dt ! time step
           integer iz ! loop index for vertical
           integer it ! loop index for time
           real time ! time from the start
           real time_pr ! time from the last print
           real time_max ! total time of the simulation
           integer ntime ! total number of timesteps
           real prtime ! frequency of output
! output variables 
           real uw(nzm+1) ! turbuelnt flux of vx
           real uwm
           real duwdz(nzm) ! vertical derivative of the turbulent flux of vx
           real vw(nzm+1) ! turbuelnt flux of vy
           real dvwdz(nzm) ! vertical derivative of the turbulent flux of vy
           real wtke(nzm+1) ! turbuelnt flux of tke
           real dwtkedz(nzm) ! vertical derivative of the turbulent flux of tke
           real wtrc(nzm+1) ! turbuelnt flux of tracer
           real dwtrcdz(nzm) ! vertical derivative of the turbulent flux of tracer
           
! flag to distinguish between alligned and staggered arrays (1: staggered, 2: alligned)
           integer iconfig
! working
           real lambdaf_x,lambdaf_y
           real tke_old(nzm) ! tke
           real dtot,wrk(8)
           real trc_an(nzm),trc_2(nzm)
           integer i

! read the input, define the mesh and obstacle paramters and intialize
           call read_input(nzm,nz,z,dz,dt,time_max,ntime,prtime,utau_x,utau_y,dpdx,dpdy, &
                           pb,ss,wx,wy,bx,by,hmean,vl,sf) 
           iconfig=1
	 if(iconfig.eq.1)write(*,*)'staggered'
	 if(iconfig.eq.2)write(*,*)'alligned'
! compute the drag coefficent and the length scales
           if(iconfig.eq.1)then
            call drag_length_stag(nzm,nz,z,wx,wy,bx,by,hmean,ceps,ck,cmu,cdragx,cdragy,dls,dlk)
           elseif(iconfig.eq.2)then
            call drag_length_alligned(nzm,nz,z,wx,wy,bx,by,hmean,ceps,ck,cmu,cdragx,cdragy,dls,dlk)
	 endif
	 write(*,*)'cdragx',cdragx(1)
! initialize vx,vy,tke - just put some value, the code should forget it since we impose a pressure gradient
           lambdaf_x=hmean*by/((wx+bx)*(wy+by))
           lambdaf_y=hmean*bx/((wx+bx)*(wy+by))
           vx=(dpdx/cdragx(1)/lambdaf_x)**.5
           vy=(dpdy/cdragy(1)/lambdaf_y)**.5
           trc=0.
           tke=0.1
           cdz=0.1
           emi=1.

           
           time=0.
           time_pr=0.
           write(*,*)'prtime',prtime
           open(unit=66,file='output')

! start the loop on time 
           do it=1,ntime
            srex_vx=0.
            srim_vx=0.
            srex_vy=0.
            srim_vy=0.
            srex_tke=0.
            srim_tke=0.
            srex_trc=0.
            srim_trc=0.
             time=time+dt
             time_pr=time_pr+dt

! first compute the diffusion coefficients

             call cdtur(nzm,nz,dz,ck,tke,dlk,cdz)

! compute the implicit and explicit terms of the equation for vx,vy and tke

             call building(nzm,nz,z,dz,pb,vl,vx,vy,cdragx,cdragy,wx,wy,bx,by,srim_vx,srim_vy,srex_tke)
             call tke_bougeault(nzm,nz,ceps,dz,vx,vy, &
                                tke,cdz,dls,sf, &                      
                                srim_tke,srex_tke)
         !    write(*,*)'srex_vx',srex_vx

! add the pressure gradient

             srex_vx=srex_vx+dpdx
             srex_vy=srex_vy+dpdy
           !  srex_vx(nz)=utau_x*utau_x

             srex_trc(1)=emi/dz/vl(1)
             srex_trc(nz)=-emi/dz
 
     !         write(*,*)'dtot=',dtot
             pr=1./0.9
             cdt=pr*cdz
! compute the vertical diffusion
            call diff(nzm,nz,1,1,dt,vx,cdz,srim_vx,srex_vx,sf,vl,dz,uw,duwdz)
            call diff(nzm,nz,1,1,dt,vy,cdz,srim_vy,srex_vy,sf,vl,dz,vw,dvwdz)
            call diff(nzm,nz,1,1,dt,tke,cdz,srim_tke,srex_tke,sf,vl,dz,wtke,dwtkedz)
            call diff(nzm,nz,1,2,dt,trc,cdt,srim_trc,srex_trc,sf,vl,dz,wtrc,dwtrcdz)
        !    write(*,*)time/3600.,trc(1)
              if(time_pr.gt.prtime)then
               time_pr=dt
               write(*,*)'time=',time,'time_pr=',time_pr
               write(66,*)'time=',time,'time_pr=',time_pr
               do iz=1,nz
                uwm=(uw(iz)+uw(iz+1))/2.
                write(66,*)iz,vx(iz)/utau_x,tke(iz)/(utau_x**2.),uwm/(utau_x**2.)
               enddo
          
              endif
           enddo
           close(66)
           trc=trc-trc(nz)
           trc_2(nz)=0.
           do iz=nz-1,1,-1
            trc_an(iz)=trc_an(iz+1)+dz/cdz(iz+1)*(emi-dwtrc(iz+1))
            trc_2(iz)=trc_2(iz+1)+dz/cdz(iz+1)*(emi)
           enddo
           if(iconfig.eq.1)open(unit=68,file='iBEPoutput_staggered_r')
           if(iconfig.eq.2)open(unit=68,file='iBEPoutput_aligned_r')
            write(68,'(a3,11(2x,a15))')'iz','vx/utau','tke/utau**2.','uwm/utau**2.','trc','wtrc','cdz','dls'
            do iz=1,nz
              uwm=(uw(iz)+uw(iz+1))/2.
              write(68,'(i3,11(2x,f15.8))')iz,vx(iz)/utau_x,tke(iz)/(utau_x**2.),uwm/(utau_x**2.),trc(iz),wtrc(iz),cdz(iz),dls(iz)
            enddo

           stop
           end   

!-------------------------------------------------------------
           subroutine read_input(nzm,nz,z,dz,dt,time_max,ntime,prtime,utau_x,utau_y,dpdx,dpdy,pb,ss,wx,wy,bx,by,hmean,vl,sf)          
           implicit none
           integer nzm
! flow variables
           real v1(nzm)   ! x component of the wind
           real v2(nzm)   ! y component of the wind
           real tke(nzm)  ! tke
           real dpdx,dpdy ! external pressure gradient 
           real utau_x,utau_y ! external pressure gradient 
! building information
           real ss(nzm+1) ! ss(iz)=probability to have a building of height z(iz)
           real pb(nzm+1) ! pb(iz)=probbaility to have a building taller or equal to z(iz)
           real wx,wy ! distance between buildings at street level in the x and y direction respectively
           real bx,by ! building dimension in the x and y direction respectively
           real hmean ! mean buiolding height
! numerical parameters
           real dz ! vertical resolution
           integer nz ! number of vertical levels
           real z(nzm+1) ! height of the faces of the numerical grid
           real dt ! time step
           integer iz,jz ! loop indexes for vertical
           integer it ! loop index for time
           integer ntime ! total number of timesteps
           real prtime ! frequency of output
           real vl(nzm) ! fraction of air in each cell
           real sf(nzm+1) ! fraction of air at the interface between cells  
           real time_max            
! Local
           character*256 line
           real hgt(nzm+1),prb(nzm+1),sstot,lambda_p
           
!
            open(unit=20,file='input_column')
10          read(20,'(A256)',END=20) line

            if (index(line,'number of points (nz)').ne.0) then
             read(20,*) nz
             if(nz.gt.nzm)then
              write(*,*)'increase nzm to ',nz
              stop
             endif
            end if

            if (index(line,'vertical resolution (dz)').ne.0) then
             read(20,*) dz
            end if

            if (index(line,'pressure forcing (utau_x, utau_y)').ne.0) then
             read(20,*) utau_x,utau_y
            end if

            if (index(line,'building dimensions (bx,by)').ne.0) then
             read(20,*) bx,by
            end if

            if (index(line,'distance between buildings (wx,wy)').ne.0) then
             read(20,*) wx,wy
            end if

            if (index(line,'building heights and probabilities (must be an integer number of dz)').ne.0) then
              jz=1
              hgt=-100.
              prb=-100.
              ss=0.
              pb=0.
14            read(20,*,err=15)hgt(jz),prb(jz)
              write(*,*)'hgt(jz),prb(jz)',hgt(jz),prb(jz)
              jz=jz+1
              goto 14
15            continue 
             
            end if
            
           if (index(line,'time step (dt)').ne.0) then
             read(20,*) dt
            end if

           if (index(line,'total time of simulation (time_max)').ne.0) then
             read(20,*) time_max
            end if

            if (index(line,'frequency of output (prtime)').ne.0) then
             read(20,*) prtime
            end if

            
      go to 10

 20   continue

      

! define the mesh

             do iz=1,nzm+1
              z(iz)=(iz-1)*dz
             enddo
! define probabilities
             do jz=1,nzm+1
               if(hgt(jz).ge.0.)then
                do iz=1,nzm+1
                 if(z(iz).eq.hgt(jz))then
                  ss(iz)=prb(jz)
                 endif
                enddo
               endif
              enddo
              do iz=1,nzm+1
               sstot=0.
               do jz=iz,nzm+1
                sstot=sstot+ss(jz)
               enddo
               pb(iz)=sstot
             enddo 
! compute fraction of air for each grid cell
             lambda_p=bx*by/((wx+bx)*(wy+by)) 
             do iz=1,nzm
              vl(iz)=1.-lambda_p*pb(iz+1)
              sf(iz)=1.-lambda_p*ss(iz)
             enddo
              sf(nzm+1)=1.
! compute mean height
             hmean=0.
             do iz=1,nz
              hmean=hmean+z(iz)*ss(iz)
             enddo
          
!
! pressure gradient
            dpdx=(utau_x**2.)/z(nz)
            dpdy=(utau_y**2.)/z(nz)
             ntime=time_max/dt

      return
      end
!-----------------------------------------------------------------
      subroutine drag_length_stag(nzm,nz,z,wx,wy,bx,by,hmean,ceps,ck,cmu,cdragx,cdragy,dls,dlk)
      implicit none
! input
      integer nzm,nz,iz
      real z(nzm+1)
      real dls(nzm),dlk(nzm)      
      real wx,wy,bx,by,hmean
      real ceps,ck,cmu
! output
      real cdragx(nzm),cdragy(nzm)
!local
      real disp,lambda_p,a1,a2,zc,d2   

! New mixing length closure (Dimitrios)
      real ellmax, zellmax, ellH, ellz
      real  kappa
      kappa = 0.41

! Here I follow the formualtion for staggered arrays based on BLM2010. To be modified for aligned

! compute the displacement height
      lambda_p=bx*by/((wx+bx)*(wy+by))
      disp=hmean*(lambda_p)**(0.13)

! Drag coefficient
      if(lambda_p.le.0.29)then
       do iz=1,nz
        cdragx(iz)=3.32*(lambda_p)**(0.47)
        cdragy(iz)=3.32*(lambda_p)**(0.47)
       enddo
      else
       do iz=1,nz
        cdragx(iz)=1.85
        cdragy(iz)=1.85
       enddo
      endif
!
      cdragx=cdragx*1.7
!

!	  write(*,*)cdragx(1),(lambda_p)**(0.47)
! length scales
      a1=2.19
      a2=1.2
      do iz=1,nz





! New mixing length closure (Dimitrios)

   ellmax   = 6.82 * lambda_p * (1.0 - lambda_p)**5
   zellmax  = 1.0 - 1.03 * lambda_p**0.2 * (1.0 - lambda_p)**1.89
   ellH     = kappa * (1.0 - lambda_p) / (1.0 + 21.32 * lambda_p)
   ellH     = kappa * (1.0 - lambda_p**0.13)

   zc       = (z(iz) + z(iz+1)) / 2.0

   if (lambda_p == 0.0) then
    ellz = kappa * zc
   else
     if (zc <= hmean) then
       if (zc/hmean <= zellmax) then
        ellz = ellmax * (zc/hmean) / zellmax
       else
         ellz = ellmax + (ellH - ellmax) * ((zc/hmean - zellmax) / (1.0 - zellmax))
       endif
     else
       ellz = ( kappa * (zc/hmean - disp/hmean))
     endif
     endif

       dlk(iz) = ellz*hmean

       ! update dls and dlk
       dls(iz) = ceps *ck  / cmu * dlk(iz)
      enddo
      
      return
      end

!----------------------------------------
! ===6================================================================72
!-----------------------------------------------------------------
      subroutine drag_length_alligned(nzm,nz,z,wx,wy,bx,by,hmean,ceps,ck,cmu,cdragx,cdragy,dls,dlk)
      implicit none
! input
      integer nzm,nz,iz
      real z(nzm+1)
      real dls(nzm),dlk(nzm)      
      real wx,wy,bx,by,hmean
      real ceps,ck,cmu
! output
      real cdragx(nzm),cdragy(nzm)
!local
      real disp,lambda_p,a1,a2,zc,d2,lambda_c,lambda_s,a,b,c,f,i,j,fs,fc,fcs      
! Here I follow the formualtion for alligned arrays 

! compute the displacement height
      lambda_p=bx*by/((wx+bx)*(wy+by))
      disp=hmean*(lambda_p)**(0.13)
! define the sheltering and channeling paramter
      lambda_c=wy/by
	  lambda_s=wx/hmean
	  write(*,*)'lambda_c=',lambda_c,'lambda_s',lambda_s
      a=0.24
	  b=1.67
	  c=2.07
	  f=0.6
	  i=1.4
	  j=4.
      fs=1.-exp(-a*(lambda_s**b))
	  fc=c/lambda_c
	  fcs=1.+f/(lambda_s**i)/(lambda_c**j)
	  
! Drag coefficient
      ! put to 0 cdragy for the moment
       do iz=1,nz
        cdragx(iz)=fc*fs*fcs
        cdragy(iz)=cdragx(iz)
       enddo
  !    write(*,*)cdragx(1),fc,fs,fcs
! length scales
      a1=2.19
      a2=1.2
      do iz=1,nz
       zc=(z(iz)+z(iz+1))/2.
       if((zc/hmean).le.1.)then
        !dls(iz)=ceps*a1*(hmean-disp)*2.
		dls(iz)=ceps*a1*(hmean-disp)
       elseif((zc/hmean).gt.1..and.(zc/hmean).le.1.5)then
        dls(iz)=ceps*a1*(zc-disp)
        
       elseif((zc/hmean).gt.1.5)then
        d2=(1.-a1/a2)*1.5*hmean+a1/a2*disp
        dls(iz)=ceps*a2*(zc-d2)
       endif
       dlk(iz)=cmu/(ceps*ck)*dls(iz)
      enddo
      
      return
      end

!----------------------------------------
! ===6================================================================72
! ===6================================================================72

      subroutine cdtur(nzm,nz,dz,ck,tke,dlk,cdz)                                  


      implicit none

! Input
! -----
      integer nzm                              ! maximum number of vertical levels
      integer nz                               ! number of vertical levels
      real dz                              ! levels size [m]
      real ck                                  ! von Karman constant
      real tke(nzm)                             ! turbulent kinetic energy
      real dlk(nzm)                             ! lenth scale

! Output
! ------
      real cdz(nzm+1)                            ! diffusion coefficient

! Local
! -----
      integer iz
      real tke_m
      real dlk_m

! ----------------------------------------------------------------------

       cdz(1)=0.

!       do iz=2,nz-1
       do iz=2,nz
        tke_m=(tke(iz-1)+tke(iz))/2.
        dlk_m=(dlk(iz-1)+dlk(iz))/2.
        cdz(iz)=ck*dlk_m*sqrt(tke_m)
       enddo
       cdz(nz+1)=cdz(nz)

       return
       end  

! ===6================================================================72

      subroutine building(nzm,nz,z,dz,pb,vl,vx,vy,cdragx,cdragy,wx,wy,bx,by,srim_vx,srim_vy,srex_tke)
      implicit none
! input      
      integer nzm
      integer nz,iz
      real z(nzm+1),dz,pb(nzm+1),ss(nzm+1),vl(nzm)
      real cdragx(nzm),cdragy(nzm),wx,wy,bx,by
      real vx(nzm),vy(nzm)
! output      
      real srim_vx(nzm),srim_vy(nzm),srex_tke(nzm)
! local
      real lfx,lfy

      do iz=1,nz
       lfx=dz*by/((wx+bx)*(wy+by))*pb(iz+1)
       lfy=dz*bx/((wx+bx)*(wy+by))*pb(iz+1)
       srim_vx(iz)=-lfx/vl(iz)/dz*cdragx(iz)*abs(vx(iz))
       srim_vy(iz)=-lfy/vl(iz)/dz*cdragy(iz)*abs(vy(iz))
! not sure it is like this for tke, but for wind orthotogonal to the face of the cube, 
! vy=0.
       srex_tke(iz)=cdragx(iz)*(lfx/vl(iz)/dz*(abs(vx(iz))**3.))+cdragy(iz)*(lfy/vl(iz)/dz*(abs(vy(iz))**3.))
      enddo

      return
      end
      
!------------------------------------------------     

      subroutine tke_bougeault (nzm,nz,ceps,dz,vx,vy, &
                                tke,cdz,dls,sf, &                       ! Input
                                srim_tke,srex_tke)                            ! Output


! ----------------------------------------------------------------------
!   Calculation of the sources (shear) and the dissipation
!    terms for the TKE equation.
! ----------------------------------------------------------------------

      implicit none

! Input
! -----
      integer nz,nzm                               ! number of vertical levels
    
      real ceps
      real dz                                  ! levels size [m]
      real vx(nzm)                              ! wind velocity
      real vy(nzm)                              ! wind velocity     
      real tke(nzm)                             ! turbulent kinetic energy
      real cdz(nzm+1)                           ! turbulent diffusion coefficient
      real dls(nzm)                             ! lrngth scale leps
      real sf(nzm+1)                            ! ?????

! Ouput
! -----
      real srim_tke(nzm)                           ! inplicit term in the tke equation
      real srex_tke(nzm)                           ! explicit term in the tke equation

! Local
! -----
      integer iz
      real sh(nzm)                              ! shear turbulent kinetic energy source
      real td(nz)                              ! dissipation

! ----------------------------------------------------------------------

      call shear_2(nzm,nz,dz,vx,vy,cdz, &                               ! Input
                  sh)                                              ! Ouput

      do iz=1,nz

       if (dls(iz).ne.0.) then
        td(iz)=-ceps*sqrt(tke(iz))/dls(iz)
       else
        td(iz)=0.
       end if ! dls
       sh(iz)=sh(iz)*sf(iz)
       srim_tke(iz)=td(iz)
       srex_tke(iz)=srex_tke(iz)+sh(iz)
      end do ! iz
      return
      end  ! 

! ===6================================================================72
! ===6================================================================72

      subroutine shear(nzm,nz,dz,vx,vy,cdz, &                               ! Input
                       sh)                                              ! Ouput



! ----------------------------------------------------------------------
!     Calculation of the shear source for the TKE equation.
! ----------------------------------------------------------------------

      implicit none

! Input
! -----
      integer nzm,nz                               ! number of vertical levels
      real dz                              ! levels size [m]
      real vx(nzm)                              ! wind velocity
      real vy(nzm)                              ! wind velocity
      real cdz(nzm+1)                           ! turbulent diffusion coefficient

! Ouput
! -----
      real sh(nzm)                              ! shear turbulent kinetic energy source

! Local
! -----
      integer iz
      real dudz1
      real dvdz1
      real dudz2
      real dvdz2
      real cdm
      real dumdz
      real cdmmin
      parameter (cdmmin=0.01)

! ----------------------------------------------------------------------

      sh(1)=0.

      do iz=2,nz-1

        dudz1=(vx(iz)-vx(iz-1))/dz
        dvdz1=(vy(iz)-vy(iz-1))/dz
        dudz2=(vx(iz+1)-vx(iz))/dz
        dvdz2=(vy(iz+1)-vy(iz))/dz

        cdm=max(0.5*(cdz(iz)+cdz(iz+1)),cdmmin)

        dumdz=0.5*((dudz1**2.+dvdz1**2)+(dudz2**2.+dvdz2**2))

        sh(iz)=cdm*dumdz
   !    write(68,*)iz,sh(iz),dudz2,vx(iz+1)

      enddo  ! iz

      sh(nz)=0.

      return ! shear
      end

! ===6================================================================72 
! ===6================================================================72

      subroutine shear_2(nzm,nz,dz,vx,vy,cdz, &                               ! Input
                       sh)                                              ! Ouput



! ----------------------------------------------------------------------
!     Calculation of the shear source for the TKE equation.
! ----------------------------------------------------------------------

      implicit none

! Input
! -----
      integer nzm,nz                               ! number of vertical levels
      real dz                              ! levels size [m]
      real vx(nzm)                              ! wind velocity
      real vy(nzm)                              ! wind velocity
      real cdz(nzm+1)                           ! turbulent diffusion coefficient

! Ouput
! -----
      real sh(nzm)                              ! shear turbulent kinetic energy source

! Local
! -----
      integer iz
      real dudz1
      real dvdz1
      real dudz2
      real dvdz2
      real cd1,cd2

     

! ----------------------------------------------------------------------

      sh(1)=0.

      do iz=2,nz-1

        dudz1=(vx(iz)-vx(iz-1))/dz
        dvdz1=(vy(iz)-vy(iz-1))/dz
        dudz2=(vx(iz+1)-vx(iz))/dz
        dvdz2=(vy(iz+1)-vy(iz))/dz

      
        sh(iz)=0.5*(cdz(iz)*(dudz1**2.+dvdz1**2)+cdz(iz+1)*(dudz2**2.+dvdz2**2))
    !   write(68,*)iz,sh(iz),dudz2,vx(iz+1)

      enddo  ! iz

      sh(nz)=0.

      return ! shear_2
      end

! ===6================================================================72 
! ===6=8===============================================================72

       subroutine diff(nzm,nz,iz1,izf,dt,co,cd,aa,bb,sf,vl,dz,fc,df)

!     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     +           Subroutine prepared by A.Martilli                     +
!     +        Ecole Polytechnique Federale de Lausanne                 +
!     +   DGR - IGE - Laboratoire de Pollution de l'Air et des Sols     +
!     +            tel.: (021)-693-61-60                                +
!     +            Email: alberto.martilli@dgr.epfl.ch                  + 
!     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!------------------------------------------------------------------------
!           Calculation of the diffusion in 1D        
!------------------------------------------------------------------------
!  - Input:
!       nz    : number of points
!       iz1   : first calculated point
!       co    : concentration of the variable of interest
!       dz    : vertical levels
!       cd    : diffusion coefficients
!       dtext : external time step
!       itest : if itest eq 1 then update co, else store in a flux array
!  - Output:
!       co :concentration of the variable of interest

!  - Internal:
!       cddz  : constant terms in the equations 
!       dt    : diffusion time step
!       nt    : number of the diffusion time steps
!       cstab : ratio of the stability condition for the time step
!---------------------------------------------------------------------

         implicit none
                
         integer nz,nzm,iz,iz1,izf
         real co(nzm),cd(nzm+1),dz,dt,dzv

         real cddz(nzm+2),fc(nzm+1),df(nzm)
         real a(nzm,3),c(nzm)
         real sf(nzm+1),vl(nzm)
         real aa(nzm),bb(nzm)
        
! Compute cddz=2*cd/dz  
        
        cddz(1)=sf(1)*cd(1)/dz
        do iz=2,nz
         cddz(iz)=2.*sf(iz)*cd(iz)/(2.*dz)
        enddo
        if(izf.gt.0)then
         cddz(nz+1)=sf(nz+1)*cd(nz+1)/dz
        else
         cddz(nz+1)=0.
        endif
         do iz=1,iz1-1
          a(iz,1)=0.
          a(iz,2)=1.
          a(iz,3)=0.
          c(iz)=co(iz)
         enddo
          
          do iz=iz1,nz-izf
           dzv=vl(iz)*dz
           a(iz,1)=-cddz(iz)*dt/dzv         
           a(iz,2)=1+dt*(cddz(iz)+cddz(iz+1))/dzv-aa(iz)*dt
           a(iz,3)=-cddz(iz+1)*dt/dzv            
           c(iz)=co(iz)+bb(iz)*dt                     
          enddo

          if(izf.eq.1)then
           dzv=vl(nz)*dz
           a(nz,1)=-cddz(nz)*dt/dzv        
           a(nz,2)=1+dt*(cddz(nz))/dzv-aa(nz)*dt
           a(nz,3)=0.            
           c(nz)=co(nz)+bb(nz)*dt   
          else
           do iz=nz-izf+1,nz
            a(iz,1)=0.
            a(iz,2)=1.
            a(iz,3)=0.
            c(iz)=co(iz)
           enddo
          endif

          call invert (nzm,nz,a,c,co)
                  
          do iz=1,iz1 
           fc(iz)=0.
          enddo
          
          do iz=iz1+1,nz 
           fc(iz)=-(cddz(iz)*(co(iz)-co(iz-1)))
          enddo
        
          do iz=1,iz1
           df(iz)=0.
          enddo
          
          do iz=iz1+1,nz-izf
           dzv=vl(iz)*dz
           if(iz.lt.nz)then
            df(iz)=+(co(iz-1)*cddz(iz)&
                 -co(iz)*(cddz(iz)+cddz(iz+1))&
                 +co(iz+1)*cddz(iz+1))/dzv
           else
            df(iz)=+(co(iz-1)*cddz(iz)&
                 -co(iz)*(cddz(iz)+cddz(iz+1)))/dzv
           endif
          enddo
          
          do iz=nz-izf,nz 
           df(iz)=0.
          enddo
                                        
       return
       end
!----------------------------------------------------------------------------

       subroutine invert(nzm,nn,a,c,x)
       implicit none
       
!ccccccccccccccccccccccccccccccc       
! Aim: INversion and resolution of a tridiagonal matrix
!          A X = C
! Input:
!  a(*,1) lower diagonal (Ai,i-1)
!  a(*,2) principal diagonal (Ai,i)
!  a(*,3) upper diagonal (Ai,i+1)
!  c      
! Output
!  x     results
!ccccccccccccccccccccccccccccccc 
       integer nzm,nn,in
       real a(nzm,3),c(nzm),x(nzm)                       
        
        do in=nn-1,1,-1                 
         c(in)=c(in)-a(in,3)*c(in+1)/a(in+1,2)
         a(in,2)=a(in,2)-a(in,3)*a(in+1,1)/a(in+1,2)
        enddo
        
        do in=2,nn        
         c(in)=c(in)-a(in,1)*c(in-1)/a(in-1,2)
        enddo
        
        do in=1,nn
         x(in)=c(in)/a(in,2)
        enddo

        return
        end
          
         

! ===6=8===============================================================72

      subroutine gaussj(a,n,b,np)

! ----------------------------------------------------------------------
! This routine solve a linear system of n equations of the form
!              A X = B
!  where  A is a matrix a(i,j)
!         B a vector and X the solution
! In output b is replaced by the solution
! ----------------------------------------------------------------------

!      include 'wind.h'

! ----------------------------------------------------------------------
! INPUT:
! ----------------------------------------------------------------------
      integer np
      real a(np,np)

! ----------------------------------------------------------------------
! OUTPUT:
! ----------------------------------------------------------------------
      real b(np)

! ----------------------------------------------------------------------
! LOCAL:
! ----------------------------------------------------------------------
      integer nm
      parameter (nm=150)

      real*8 big,dum
      integer i,icol,irow
      integer j,k,l,ll,n
      integer ii,jj
      integer ipiv(nm)
      real*8 pivinv

! ----------------------------------------------------------------------
! END VARIABLES DEFINITIONS
! ----------------------------------------------------------------------


      do j=1,n
         ipiv(j)=0.
      enddo

      do i=1,n
         big=0.
         do j=1,n
            if(ipiv(j).ne.1)then
               do k=1,n
                  if(ipiv(k).eq.0)then
                     if(abs(a(j,k)).ge.big)then
                        big=abs(a(j,k))
                        irow=j
                        icol=k
                     endif
                  elseif(ipiv(k).gt.1)then
                     !pause'singular matrix in gaussj'
                 print *, 'ERROR: Singular matrix in gaussj'
                  stop 1

                  endif
               enddo
            endif
         enddo

         ipiv(icol)=ipiv(icol)+1

         if(irow.ne.icol)then
            do l=1,n
               dum=a(irow,l)
               a(irow,l)=a(icol,l)
               a(icol,l)=dum
            enddo

            dum=b(irow)
            b(irow)=b(icol)
            b(icol)=dum

         endif

         !if(a(icol,icol).eq.0)pause 'singular matrix in gaussj'
         if(a(icol,icol).eq.0) then
         print *, 'ERROR: Singular matrix in gaussj'
         stop 1
         endif

         pivinv=1./a(icol,icol)
         a(icol,icol)=1

         do l=1,n
            a(icol,l)=a(icol,l)*pivinv
         enddo

         b(icol)=b(icol)*pivinv

         do ll=1,n
            if(ll.ne.icol)then
               dum=a(ll,icol)
               a(ll,icol)=0.
               do l=1,n
                  a(ll,l)=a(ll,l)-a(icol,l)*dum
               enddo

               b(ll)=b(ll)-b(icol)*dum

            endif
         enddo


      enddo

      return
      end !subroutine gaussj

! ===6=8===============================================================72




      
         

