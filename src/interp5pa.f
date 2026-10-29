      program interp5pa
c=====================================================================
c  Rasterise the five HIST-ETAS parameter fields onto a regular grid.
c
c  hist-etas5pa represents each varying parameter as
c
c        theta(x,y) = theta_0 * exp( f_theta(x,y) )
c
c  with f_theta piecewise linear on the Delaunay tessellation held in
c  delone2.out and its vertex coefficients in hist-etas5pa.upda.  The
c  package plots those coefficients as coloured dots at the vertices
c  (Section11.5files/hist-etas5pa.R); this program evaluates the field
c  itself on a lon/lat grid so it can be drawn as a continuous map.
c
c  Interpolation rule, taken from the package rather than invented here:
c    mu   -- linear in the VALUE.  hist-etas5pa.f integrates the
c            background over a triangle as b1(idk)/3*adt(id)*(tz-tstart),
c            i.e. the mean of the three vertex rates times the area,
c            which is the integral of the value-linear interpolant.  The
c            conditional-intensity program interpolates mu the same way.
c    K, alpha, p, q -- linear in the LOG, i.e. in f_theta.  That is what
c            workpara5pa.f does when it needs these off-vertex.
c
c  alpha is reported as a4*exp(f), the full Utsu productivity exponent,
c  matching hist-etas5pa.R.  The kernel itself uses a4/2 because the
c  scaling enters the spatial variance: exp(2*(a4/2)*(M-M0)).
c
c  Inputs   node.conf                 lon/lat window and pixel counts
c           hist-etas5pa-forecast.conf   region origin (xmin, ymin)
c           delone2.out               vertices + triangles
c           hist-etas5pa.upda         base parameters + vertex field
c  Output   interp5pa.pixel   lon lat mu K alpha p q   (-999 = outside
c                             the tessellation's convex hull)
c=====================================================================
      implicit real * 8 (a-h,o-z)
      parameter(npa=8000, mdl=16000)
      dimension xx(npa),yy(npa),b(npa*5)
      dimension flog(npa,5),ref(5)
      dimension ldt(mdl,3)
      dimension tbxmin(mdl),tbxmax(mdl),tbymin(mdl),tbymax(mdl)
      dimension xd(3),yd(3),val(5)
      det(x1,x2,x3,y1,y2,y3)=x1*y2+x2*y3+x3*y1-x1*y3-x2*y1-x3*y2
      eps=1.d-10

c---- grid window -----------------------------------------------------
      open(1,file='node.conf')
      read(1,*) glon1,glon2,glat1,glat2
      read(1,*) mx,my
      close(1)

c---- region origin ---------------------------------------------------
      open(1,file='hist-etas5pa-forecast.conf')
      read(1,*) tx,ty,tz
      read(1,*) xxmin,yymin
      close(1)

c---- tessellation ----------------------------------------------------
      open(1,file='delone2.out')
      read(1,*) nn,npext,ndelo,tx,ty
      if(npext.gt.npa) stop 'interp5pa: npext exceeds npa'
      if(ndelo.gt.mdl) stop 'interp5pa: ndelo exceeds mdl'
      do 10 i=1,npext
        read(1,*) idum,xx(i),yy(i),am,ot,dd,ee,ff
   10 continue
      do 20 i=1,ndelo
        read(1,*) idum,(ldt(i,j),j=1,3),adum
   20 continue
      close(1)

c---- fitted fields ---------------------------------------------------
      open(1,file='hist-etas5pa.upda')
      read(1,*) w1,w2,w3,w4,w5,abic,m,w0
      read(1,*) a1,a2,a3,a4,a5,a6,a7
      read(1,*) (b(i),i=1,m)
      close(1)
      m5=m/5
      if(m5.ne.npext) then
        write(6,*) 'interp5pa: vertex count mismatch',m5,npext
        stop
      endif
      ref(1)=a1
      ref(2)=a2
      ref(3)=a4
      ref(4)=a5
      ref(5)=a7
      do 30 i=1,m5
        do 30 k=1,5
          flog(i,k)=b(i+(k-1)*m5)
   30 continue
      write(6,*) 'npext,ndelo,abic',npext,ndelo,abic
      write(6,*) 'mu0,K0,c,alpha,p,d,q'
      write(6,*) a1,a2,a3,a4,a5,a6,a7

c---- triangle bounding boxes ----------------------------------------
      do 40 ipt=1,ndelo
        bx1= 1.d69
        bx2=-1.d69
        by1= 1.d69
        by2=-1.d69
        do 35 k=1,3
          xk=xx(ldt(ipt,k))
          yk=yy(ldt(ipt,k))
          if(xk.lt.bx1) bx1=xk
          if(xk.gt.bx2) bx2=xk
          if(yk.lt.by1) by1=yk
          if(yk.gt.by2) by2=yk
   35   continue
        tbxmin(ipt)=bx1
        tbxmax(ipt)=bx2
        tbymin(ipt)=by1
        tbymax(ipt)=by2
   40 continue

c---- rasterise -------------------------------------------------------
      open(9,file='interp5pa.pixel')
      nin=0
      nout=0
      do 100 i=1,mx
      do 100 j=1,my
        alon=glon1+(glon2-glon1)*(i-0.5d0)/mx
        alat=glat1+(glat2-glat1)*(j-0.5d0)/my
        x=alon-xxmin
        y=alat-yymin
        do 110 ipt=1,ndelo
          if(x.lt.tbxmin(ipt)) go to 110
          if(x.gt.tbxmax(ipt)) go to 110
          if(y.lt.tbymin(ipt)) go to 110
          if(y.gt.tbymax(ipt)) go to 110
          do 115 k=1,3
            xd(k)=xx(ldt(ipt,k))
            yd(k)=yy(ldt(ipt,k))
  115     continue
          denom=det(xd(1),xd(2),xd(3),yd(1),yd(2),yd(3))
          if(denom.eq.0.0d0) go to 110
          ax1=det(x    ,xd(2),xd(3),y    ,yd(2),yd(3))/denom
          ax2=det(xd(1),x    ,xd(3),yd(1),y    ,yd(3))/denom
          ax3=det(xd(1),xd(2),x    ,yd(1),yd(2),y    )/denom
          if(ax1.lt.eps) go to 110
          if(ax2.lt.eps) go to 110
          if(ax3.lt.eps) go to 110
          j1=ldt(ipt,1)
          j2=ldt(ipt,2)
          j3=ldt(ipt,3)
c         mu: linear in the value (matches the likelihood's triangle
c         integral and the conditional-intensity program)
          val(1)=ax1*ref(1)*exp(flog(j1,1))
     &          +ax2*ref(1)*exp(flog(j2,1))
     &          +ax3*ref(1)*exp(flog(j3,1))
c         K, alpha, p, q: linear in the log (matches workpara5pa.f)
          do 120 k=2,5
            val(k)=ref(k)*exp(ax1*flog(j1,k)
     &                       +ax2*flog(j2,k)
     &                       +ax3*flog(j3,k))
  120     continue
          nin=nin+1
          go to 130
  110   continue
        do 125 k=1,5
          val(k)=-999.d0
  125   continue
        nout=nout+1
  130   continue
        write(9,1) alon,alat,(val(k),k=1,5)
    1   format(2f10.4,5e15.6)
  100 continue
      close(9)
      write(6,*) 'pixels inside/outside hull:',nin,nout
      stop
      end
