!
#include <defs.h>      
      module ringfilter_module
!
! Coded by Tong Dang, May, 2017
! This program is to build a new ring filter for the future high-resolution 
! TIE-GCM. The filter procedures are as follows:
! 1. Divde the longitude sessions into several trunks
! 2. Calculate the average values of each trunck
! 3. Do the restruction based on the mean value of each trunk
!
      use params_module,only: nlon,nlonp4,nlat
      implicit none
!
      contains
!-----------------------------------------------------------------------
      subroutine ringfilter(f,lev0,lev1,lat)
!
! Ringfilter for the second order of FFT
! keep first and second order of fourier series, and filter orders
! coded by Dang, 2017
! Args:
      real,intent(inout) :: f(nlonp4,lev0:lev1)
      integer,intent(in) :: lat,lev0,lev1
!
! Local:
      real :: fx(nlon),average(200),f_out(nlon),
     |  w(nlon),wm1(nlon),a0,a1,b1,theta,dtheta
     |  ,fL,fR,fm2,fm1,ff,fp1,fp2,a,b,c
#if (NLAT==36 && NLON==72)
      integer,parameter :: nlat_filter = 4
#elif (NLAT==72 && NLON==144)
      integer,parameter :: nlat_filter = 6
#elif (NLAT==144 && NLON==288) 
      integer,parameter :: nlat_filter = 18
#elif (NLAT==288 && NLON==576)
      integer,parameter :: nlat_filter = 40
#else
      UNKNOWN NLAT,NLON      ! compilation will stop here if unknown res
#endif
      integer :: chunk_array(nlat_filter),i,j,jj,k,
     |  points,n,chunk,j_index,m
!
#if (NLAT==36 && NLON==72)
      chunk_array=(/9,9,18,36/) 
#elif (NLAT==72 && NLON==144)
      chunk_array=(/9,18,36,36,72,72/) ! 2.5 deg 8 intotal
#elif (NLAT==144 && NLON==288)
      chunk_array=(/9,9,18,18,36,36,36,36,72,72,144,144,144,144,
     | 144,144,144,144/) ! 1.25 deg 18 intotal
#elif (NLAT==288 && NLON==576)
      chunk_array=(/9,9,9,9,18,18,18,18,36,36,36,36,36,36,36,36,
     | 72,72,72,72,72,72,72,72,144,144,144,144,144,144,144,144,
     | 288,288,288,288,288,288,288,288/) ! 0.625 deg 40 intotal
#else
      UNKNOWN NLAT,NLON      ! compilation will stop here if unknown res
#endif        
        
      if(lat .LE. nlat_filter .OR. lat .GE. (nlat-nlat_filter+1) ) then

      dtheta=2.*3.14159/real(nlon)

      do k=lev0,lev1

! Load field data into w
! Fourier expansion: f(x)=a0+a1*cos(x)+b1*sin(x)+others
        a1=0.
        b1=0.
        do i=1,nlon
           w(i) = f(i+2,k)
           theta=dtheta*i
           a1=a1+w(i)*cos(theta)
           b1=b1+w(i)*sin(theta)
        enddo
        a1=a1*2./real(nlon)
        b1=b1*2./real(nlon)
        a0=sum(w)/real(nlon)

! Chunk numbers in this latitude
        if(lat .LE. nlat_filter) chunk=chunk_array(lat)
        if(lat .GE. (nlat-nlat_filter+1)) chunk=chunk_array(nlat-lat+1)

! w(i)=wm1(i)+fx(i), then filter fx(i)
        do i=1,nlon
          theta=dtheta*i
          wm1(i)=a0+a1*cos(theta)+b1*sin(theta)
          fx(i)=w(i)-wm1(i)
          fx(i)=w(i)
        enddo
       
! Start the ring average filtering

! Grid points in each chunk
        points=nlon/chunk
        n=points

! Calculate the average value in each chunk
        do i=1,chunk    ! i is the chunk number in each ring
          average(i)=sum(fx((i-1)*points+1:i*points))/real(points)
        enddo
        
! Then do the linear interpolation between each fL, fR
        do i=1,chunk  ! i is the chunk number in each ring

! Calculate f,fL,fR 
          if(i .EQ. 1) then

            fm2=average(chunk-1)
            fm1=average(chunk)
             ff=average(i)
            fp1=average(i+1)
            fp2=average(i+2)

          else if(i .EQ. 2) then

            fm2=average(chunk)
            fm1=average(i-1)
             ff=average(i)
            fp1=average(i+1)
            fp2=average(i+2)

          else if(i .EQ. chunk-1) then

            fm2=average(i-2)
            fm1=average(i-1)
             ff=average(i)
            fp1=average(i+1)
            fp2=average(1)

          else if(i .EQ. chunk) then

            fm2=average(i-2)
            fm1=average(i-1)
             ff=average(i)
            fp1=average(1)
            fp2=average(2)

          else

            fm2=average(i-2)
            fm1=average(i-1)
             ff=average(i)
            fp1=average(i+1)
            fp2=average(i+2)

          endif

          fL = (-fm2+7.*fm1+7.*ff-fp1)/12.
          fR = (-fm1+7.*ff+7.*fp1-fp2)/12.

          a = 3.*(fL + fR - 2.*ff)
          b = 2.*(3.*ff - fR - 2.*fL)
          c = fL

! Calculate the filtered data at j_index
          do m=1,n
            j_index=m+(i-1)*points
!            f_out(j_index)=fL*(2.*n-2.*m+1.)/(2.*n)+
!     |       fR*(2.*m-1.)/(2.*n)
            f_out(j_index)=(a/3.0)*(3*m*m-3*m+1)/(n*n) 
     |       + 0.5*b*(2*m-1)/n + c
          enddo

        enddo ! i=1,chunk

        fx(:)=f_out(:)

! Save filtered field:
        do i=1,nlon
          f(i+2,k) = fx(i) + wm1(i)
          f(i+2,k) = fx(i)
        enddo ! i=1,nlon

       enddo ! k=lev0,lev1

      endif ! lat
      end subroutine ringfilter
!-----------------------------------------------------------------------
      subroutine ringfilter_uv(f,lev0,lev1,lat)
!
! Ringfilter for the second order of FFT
! keep first and second order of fourier series, and filter orders
! coded by Dang, 2017
! Args:
      real,intent(inout) :: f(nlonp4,lev0:lev1)
      integer,intent(in) :: lat,lev0,lev1
!
! Local:
      real :: fx(nlon),average(200),f_out(nlon),
     |  w(nlon),wm1(nlon),a0,a1,b1,a2,b2,a3,b3,theta,dtheta
     |  ,fL,fR,fm2,fm1,ff,fp1,fp2,a,b,c
!
#if (NLAT==36 && NLON==72)
      integer,parameter :: nlat_filter = 4
#elif (NLAT==72 && NLON==144)
      integer,parameter :: nlat_filter = 6
#elif (NLAT==144 && NLON==288) 
      integer,parameter :: nlat_filter = 18
#elif (NLAT==288 && NLON==576)
      integer,parameter :: nlat_filter = 40
#else
      UNKNOWN NLAT,NLON      ! compilation will stop here if unknown res
#endif
      integer :: chunk_array(nlat_filter),i,j,jj,k,
     |  points,n,chunk,j_index,m
!
#if (NLAT==36 && NLON==72)
      chunk_array=(/9,9,18,36/)
#elif (NLAT==72 && NLON==144)
      chunk_array=(/9,18,36,36,72,72/) ! 2.5 deg 8 intotal
#elif (NLAT==144 && NLON==288)
      chunk_array=(/9,9,18,18,36,36,36,36,72,72,144,144,144,144,
     | 144,144,144,144/) ! 1.25 deg 18 intotal
#elif (NLAT==288 && NLON==576)
      chunk_array=(/9,9,9,9,18,18,18,18,36,36,36,36,36,36,36,36,
     | 72,72,72,72,72,72,72,72,144,144,144,144,144,144,144,144,
     | 288,288,288,288,288,288,288,288/) ! 0.625 deg 40 intotal
#else
      UNKNOWN NLAT,NLON      ! compilation will stop here if unknown res
#endif        

      if(lat .LE. nlat_filter .OR. lat .GE. (nlat-nlat_filter+1) ) then

      dtheta=2.*3.14159/real(nlon)

      do k=lev0,lev1

! Load field data into w
! Fourier expansion: f(x)=a0+a1*cos(x)+b1*sin(x)+a2*cos(2x)+b2*sin(2x)+
!                        a3*cos(3x)+b3*sin(3x)
        a1=0.
        b1=0.
        a2=0.
        b2=0.
        a3=0.
        b3=0.
        do i=1,nlon
           w(i) = f(i+2,k)
           theta=dtheta*i
           a1=a1+w(i)*cos(theta)
           b1=b1+w(i)*sin(theta)
           a2=a2+w(i)*cos(2*theta)
           b2=b2+w(i)*sin(2*theta)
           a3=a3+w(i)*cos(3*theta)
           b3=b3+w(i)*sin(3*theta)
        enddo
        a1=a1*2./real(nlon)
        b1=b1*2./real(nlon)
        a2=a2*2./real(nlon)
        b2=b2*2./real(nlon)
        a3=a3*2./real(nlon)
        b3=b3*2./real(nlon)
        a0=sum(w)/real(nlon)

! Chunk numbers in this latitude
        if(lat .LE. nlat_filter) chunk=chunk_array(lat)
        if(lat .GE. (nlat-nlat_filter+1)) chunk=chunk_array(nlat-lat+1)

! w(i)=wm1(i)+fx(i), then filter fx(i)
        do i=1,nlon
          theta=dtheta*i
          wm1(i)=a0+a1*cos(theta)+b1*sin(theta)+a2*cos(2*theta)+
     |     b2*sin(2*theta)+a3*cos(3*theta)+b3*sin(3*theta)
          fx(i)=w(i)-wm1(i)
        enddo
       
! Start the ring average filtering

! Grid points in each chunk
        points=nlon/chunk
        n=points

! Calculate the average value in each chunk
        do i=1,chunk    ! i is the chunk number in each ring
          average(i)=sum(fx((i-1)*points+1:i*points))/real(points)
        enddo
       
! Then do the linear interpolation between each fL, fR
        do i=1,chunk-1  ! i is the chunk number in each ring
 
! Calculate f,fL,fR 
          if(i .EQ. 1) then

            fm2=average(chunk-1)
            fm1=average(chunk)
             ff=average(i)
            fp1=average(i+1)
            fp2=average(i+2)

          else if(i .EQ. 2) then

            fm2=average(chunk)
            fm1=average(i-1)
             ff=average(i)
            fp1=average(i+1)
            fp2=average(i+2)

          else if(i .EQ. chunk-1) then

            fm2=average(i-2)
            fm1=average(i-1)
             ff=average(i)
            fp1=average(i+1)
            fp2=average(1)

          else if(i .EQ. chunk) then

            fm2=average(i-2)
            fm1=average(i-1)
             ff=average(i)
            fp1=average(1)
            fp2=average(2)

          else

            fm2=average(i-2)
            fm1=average(i-1)
             ff=average(i)
            fp1=average(i+1)
            fp2=average(i+2)

          endif

          fL = (-fm2+7.*fm1+7.*ff-fp1)/12.
          fR = (-fm1+7.*ff+7.*fp1-fp2)/12.

          a = 3.*(fL + fR - 2.*ff)
          b = 2.*(3.*ff - fR - 2.*fL)
          c = fL

! Calculate the filtered data at j_index
          do m=1,n
            j_index=m+(i-1)*points
!            f_out(j_index)=fL*(2.*n-2.*m+1.)/(2.*n)+
!     |       fR*(2.*m-1.)/(2.*n)
            f_out(j_index)=(a/3.0)*(3*m*m-3*m+1)/(n*n) 
     |       + 0.5*b*(2*m-1)/n + c
          enddo

        enddo ! i=1,chunk

        fx(:)=f_out(:)

! Save filtered field:
        do i=1,nlon
          f(i+2,k) = fx(i) + wm1(i)
        enddo ! i=1,nlon

       enddo ! k=lev0,lev1

      endif ! lat
      end subroutine ringfilter_uv
!-----------------------------------------------------------------------
      subroutine linearfilter(f,lev0,lev1,lat)
!
! Ringfilter for the second order of FFT
! keep first and second order of fourier series, and filter orders
! coded by Dang, 2017
! Args:
      real,intent(inout) :: f(nlonp4,lev0:lev1)
      integer,intent(in) :: lat,lev0,lev1
!
! Local:
      real :: fx(nlon),average(200),f_out(nlon),
     |  w(nlon),wm1(nlon),a0,a1,b1,theta,dtheta

#if (NLAT==36 && NLON==72)
      integer,parameter :: nlat_filter = 4
#elif (NLAT==72 && NLON==144)
      integer,parameter :: nlat_filter = 6
#elif (NLAT==144 && NLON==288) 
      integer,parameter :: nlat_filter = 18
#elif (NLAT==288 && NLON==576)
      integer,parameter :: nlat_filter = 40
#else
      UNKNOWN NLAT,NLON      ! compilation will stop here if unknown res
#endif
      integer :: chunk_array(nlat_filter),i,j,jj,k,
     |  points,n,chunk,j_index,m
!
#if (NLAT==36 && NLON==72)
      chunk_array=(/9,9,18,36/)
#elif (NLAT==72 && NLON==144)
      chunk_array=(/9,18,36,36,72,72/) ! 2.5 deg 8 intotal
#elif (NLAT==144 && NLON==288)
      chunk_array=(/9,9,18,18,36,36,36,36,72,72,144,144,144,144,
     | 144,144,144,144/) ! 1.25 deg 18 intotal
#elif (NLAT==288 && NLON==576)
      chunk_array=(/9,9,9,9,18,18,18,18,36,36,36,36,36,36,36,36,
     | 72,72,72,72,72,72,72,72,144,144,144,144,144,144,144,144,
     | 288,288,288,288,288,288,288,288/) ! 0.625 deg 40 intotal
#else
      UNKNOWN NLAT,NLON      ! compilation will stop here if unknown res
#endif        

      if(lat .LE. nlat_filter .OR. lat .GE. (nlat-nlat_filter+1) ) then

      dtheta=2.*3.14159/real(nlon)

      do k=lev0,lev1

! Load field data into w
! Fourier expansion: f(x)=a0+a1*cos(x)+b1*sin(x)+others
        a1=0.
        b1=0.
        do i=1,nlon
           w(i) = f(i+2,k)
           theta=dtheta*i
           a1=a1+w(i)*cos(theta)
           b1=b1+w(i)*sin(theta)
        enddo
        a1=a1*2./real(nlon)
        b1=b1*2./real(nlon)
        a0=sum(w)/real(nlon)

! Chunk numbers in this latitude
        if(lat .LE. nlat_filter) chunk=chunk_array(lat)
        if(lat .GE. (nlat-nlat_filter+1)) chunk=chunk_array(nlat-lat+1)

! w(i)=wm1(i)+fx(i), then filter fx(i)
        do i=1,nlon
          theta=dtheta*i
          wm1(i)=a0+a1*cos(theta)+b1*sin(theta)
          fx(i)=w(i)-wm1(i)
        enddo

! Start the ring average filtering
! Grid points in each chunk
        points=nlon/chunk
        n=points

! Calculate the average value in each chunk
        do i=1,chunk    ! i is the chunk number in each ring
          average(i)=sum(fx((i-1)*points+1:i*points))/real(points)
        enddo

! Then do the linear interpolation between each average point
        do i=1,chunk-1  ! i is the chunk number in each ring

! Calculate first n/2 and last n/2  
          if(i .EQ. 1) then

          do j=1,points/2

            j_index=nlon-points/2+j
            f_out(j_index)=average(chunk)*(2.*n-2.*j+1.)/(2.*n)+
     |       average(i)*(2.*j-1.)/(2.*n)

            j_index=j
            jj=j+n/2.
            f_out(j_index)=average(chunk)*(2.*n-2.*jj+1.)/(2.*n)+
     |         average(i)*(2.*jj-1.)/(2.*n)
            enddo

           endif

! Calculate the filtered data at grid of j_index
          do j=1,points    ! j is the number index in each chunk

           j_index=(i-1)*points+points/2+j
           f_out(j_index)=average(i)*(2.*n-2.*j+1.)/(2.*n)+
     |       average(i+1)*(2.*j-1.)/(2.*n)
          enddo

        enddo ! i=1,chunk

        fx(:)=f_out(:)

! Save filtered field:
        do i=1,nlon
          f(i+2,k) = fx(i) + wm1(i)
        enddo ! i=1,nlon

       enddo ! k=lev0,lev1

      endif ! lat
      end subroutine linearfilter
!-----------------------------------------------------------------------
      subroutine linearfilter_uv(f,lev0,lev1,lat)
!
! Ringfilter for the second order of FFT
! keep first and second order of fourier series, and filter orders
! coded by Dang, 2017
! Args:
      real,intent(inout) :: f(nlonp4,lev0:lev1)
      integer,intent(in) :: lat,lev0,lev1
!
! Local:
      real :: fx(nlon),average(200),f_out(nlon),
     |  w(nlon),wm1(nlon),a0,a1,b1,a2,b2,a3,b3,theta,dtheta

#if (NLAT==36 && NLON==72)
      integer,parameter :: nlat_filter = 4
#elif (NLAT==72 && NLON==144)
      integer,parameter :: nlat_filter = 6
#elif (NLAT==144 && NLON==288) 
      integer,parameter :: nlat_filter = 18
#elif (NLAT==288 && NLON==576)
      integer,parameter :: nlat_filter = 40
#else
      UNKNOWN NLAT,NLON      ! compilation will stop here if unknown res
#endif
      integer :: chunk_array(nlat_filter),i,j,jj,k,
     |  points,n,chunk,j_index,m
!
#if (NLAT==36 && NLON==72)
      chunk_array=(/9,9,18,36/)
#elif (NLAT==72 && NLON==144)
      chunk_array=(/9,18,36,36,72,72/) ! 2.5 deg 8 intotal
#elif (NLAT==144 && NLON==288)
      chunk_array=(/9,9,18,18,36,36,36,36,72,72,144,144,144,144,
     | 144,144,144,144/) ! 1.25 deg 18 intotal
#elif (NLAT==288 && NLON==576)
      chunk_array=(/9,9,9,9,18,18,18,18,36,36,36,36,36,36,36,36,
     | 72,72,72,72,72,72,72,72,144,144,144,144,144,144,144,144,
     | 288,288,288,288,288,288,288,288/) ! 0.625 deg 40 intotal
#else
      UNKNOWN NLAT,NLON      ! compilation will stop here if unknown res
#endif        

      if(lat .LE. nlat_filter .OR. lat .GE. (nlat-nlat_filter+1) ) then

      dtheta=2.*3.14159/real(nlon)

      do k=lev0,lev1

! Load field data into w
! Fourier expansion: f(x)=a0+a1*cos(x)+b1*sin(x)+a2*cos(2x)+b2*sin(2x)+
!                        a3*cos(3x)+b3*sin(3x)
        a1=0.
        b1=0.
        a2=0.
        b2=0.
        a3=0.
        b3=0.
        do i=1,nlon
           w(i) = f(i+2,k)
           theta=dtheta*i
           a1=a1+w(i)*cos(theta)
           b1=b1+w(i)*sin(theta)
           a2=a2+w(i)*cos(2*theta)
           b2=b2+w(i)*sin(2*theta)
           a3=a3+w(i)*cos(3*theta)
           b3=b3+w(i)*sin(3*theta)
        enddo
        a1=a1*2./real(nlon)
        b1=b1*2./real(nlon)
        a2=a2*2./real(nlon)
        b2=b2*2./real(nlon)
        a3=a3*2./real(nlon)
        b3=b3*2./real(nlon)
        a0=sum(w)/real(nlon)

! Chunk numbers in this latitude
        if(lat .LE. nlat_filter) chunk=chunk_array(lat)
        if(lat .GE. (nlat-nlat_filter+1)) chunk=chunk_array(nlat-lat+1)

! w(i)=wm1(i)+fx(i), then filter fx(i)
        do i=1,nlon
          theta=dtheta*i
          wm1(i)=a0+a1*cos(theta)+b1*sin(theta)+a2*cos(2*theta)+
     |     b2*sin(2*theta)+a3*cos(3*theta)+b3*sin(3*theta)
          fx(i)=w(i)-wm1(i)
        enddo

! Start the ring average filtering

! Grid points in each chunk
        points=nlon/chunk
        n=points

! Calculate the average value in each chunk
        do i=1,chunk    ! i is the chunk number in each ring
          average(i)=sum(fx((i-1)*points+1:i*points))/real(points)
        enddo

! Then do the linear interpolation between each average point
        do i=1,chunk-1  ! i is the chunk number in each ring

! Calculate first n/2 and last n/2  
          if(i .EQ. 1) then

          do j=1,points/2

            j_index=nlon-points/2+j
            f_out(j_index)=average(chunk)*(2.*n-2.*j+1.)/(2.*n)+
     |       average(i)*(2.*j-1.)/(2.*n)

            j_index=j
            jj=j+n/2.
            f_out(j_index)=average(chunk)*(2.*n-2.*jj+1.)/(2.*n)+
     |         average(i)*(2.*jj-1.)/(2.*n)
            enddo

           endif

! Calculate the filtered data at grid of j_index
          do j=1,points    ! j is the number index in each chunk

           j_index=(i-1)*points+points/2+j
           f_out(j_index)=average(i)*(2.*n-2.*j+1.)/(2.*n)+
     |       average(i+1)*(2.*j-1.)/(2.*n)
          enddo

        enddo ! i=1,chunk

        fx(:)=f_out(:)

! Save filtered field:
        do i=1,nlon
          f(i+2,k) = fx(i) + wm1(i)
        enddo ! i=1,nlon

       enddo ! k=lev0,lev1

      endif ! lat
      end subroutine linearfilter_uv
!-----------------------------------------------------------------------
      end module ringfilter_module
