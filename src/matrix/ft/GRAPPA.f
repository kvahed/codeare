!     
!     @brief            Compute GRAPPA source and target weights
!     
!     @param[in]   aln  Measured ACS lines
!     @param[in]   asz  Dimensions of ACS data 
!     @param[in]   af   Acceleration factor
!     @param[in]   ksz  Kernel size
!     @param[in]   msz  Measurement size
!     @param[in]   dim  Scan dimension (2D/3D)
!     @param[out]  w    Weights 
!
     
      SUBROUTINE d_src_trg_mat (aln, asz, msz, dim, af, s, ssz, t, tsz) 
      integer, dimension (3), intent (in) :: asz
      integer, dimension (2), intent (in) :: ssz, tsz
      integer, dimension (4), intent (in) :: msz
      integer, intent (in) :: af, dim
      complex*16, dimension (asz(1), asz(2), asz(3)), intent (in) :: aln
      complex*16, dimension (ssz(1), ssz(2)) :: s
      complex*16, dimension (tsz(1), tsz(2)) :: t
      integer :: nx, ny, nz, nc, sx, sy, dx, dy, nxacs, nyacs, nzacs
      integer :: i, x, y, nss
      
      nx = msz (1)
      ny = msz (2)
      if      (dim == 3) then
         nz = msz (3)
         nc = msz (4)
      else if (dim == 2) then
         nz = 1
         nc = msz (3)
      end if
      
      nxacs = asz(1)
      nyacs = asz(2)
      nzacs = asz(3)
      
      sx = ksz(1)
      sy = ksz(2)
      
      dx = floor(real(sx)/2.0)
      dy = (sy/2-1)*af 
      
      i = 0
      
      do x = dx+1,nxacs-dx
         do y = 1,nyacs-(sy-1)*af
            i = i+1
            s (:,i) = reshape(aln(:,y:af:y+(sy-1)*af,x-dx:x+dx),(
     $           /integer(nc*srcy*srcx)/))
            t (:,i) = reshape(aln(:,y+1+dy:y+dy+af-1,x),(/integer(nc*(af
     $           -1))/))
         end do
      end do

      end SUBROUTINE d_src_trg_mat



      SUBROUTINE s_src_trg_mat (aln, asz, msz, dim, af, s, ssz, t, tsz) 
      integer, dimension (3), intent (in) :: asz
      integer, dimension (2), intent (in) :: ssz, tsz
      integer, dimension (4), intent (in) :: msz
      integer, intent (in) :: af, dim
      complex*8 , dimension (asz(1), asz(2), asz(3)) :: aln
      complex*8 , dimension (ssz(1), ssz(2)) :: s
      complex*8 , dimension (tsz(1), tsz(2)) :: t
      integer :: nx, ny, nz, nc, sx, sy, dx, dy, nxacs, nyacs, nzacs
      integer :: i, x, y, nss
      
      nx = msz (1)
      ny = msz (2)
      if      (dim == 3) then
         nz = msz (3)
         nc = msz (4)
      else if (dim == 2) then
         nz = 1
         nc = msz (3)
      end if

      nxacs = asz(1)
      nyacs = asz(2)
      nzacs = asz(3)
      
      sx = ksz(1)
      sy = ksz(2)
      
      dx = floor(real(sx)/2.0)
      dy = (sy/2-1)*af 
      
      i = 0
      
      do x = dx+1,nxacs-dx
         do y = 1,nyacs-(sy-1)*af
            i = i+1
            s (:,i) = reshape(aln(:,y:af:y+(sy-1)*af,x-dx:x+dx),(
     $           /integer(nc*srcy*srcx)/))
            t (:,i) = reshape(aln(:,y+1+dy:y+dy+af-1,x),(/integer(nc*(af
     $           -1))/))
         end do
      end do

      end SUBROUTINE s_src_trg_mat


      SUBROUTINE d_apply_weights (n, d, af, meas, recon)
      integer, dimension (4), intent (in) :: n
      integer, dimension (2), intent (in) :: d
      integer, intent (in) :: af
c$$$      complex*16, dimension (nx,ny), intent (in) :: meas
c$$$      complex*16, 
      
      
      end SUBROUTINE d_apply_weights

c$$$      tic; fprintf('Applyig weights ')
c$$$      
c$$$      sigrecon = zeros(nc,ny+2*dy+1,nx+2*dx);
c$$$     $     % prepare matrix for convolution     
c$$$      sigrecon(:,dy+1:af:ny+dy,dx+1:end-dx) = sig;
c$$$     $     % Write undersampled data into zero-padded matrix 
c$$$      
c$$$      for xind = dx+1:nx+dx, 
c$$$      for yind= 1:af:ny,
c$$$      src=reshape(sigrecon(:,yind:af:yind + (srcy-1)*af,xind-dx:xind+dx)
c$$$     $     ,nc*srcy*srcx,1);
c$$$      sigrecon(:,yind+dy+1:yind+dy+af-1,xind)=reshape(ws*src,[nc (af
c$$$     $     -1)]);           %Apply weights to source points
c$$$      end
c$$$      end

c$$$      sigrecon = sigrecon(:,dy+1:ny+dy,dx+1:nx+dx);
c$$$     $     %Crop out the good data.

