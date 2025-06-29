c
c User subroutine VUSDFLD for user-defined fields
c
      subroutine vusdfld(
c Read only -
     *   nblock, nstatev, nfieldv, nprops, ndir, nshr, 
     *   jElemUid, kIntPt, kLayer, kSecPt, 
     *   stepTime, totalTime, dt, cmname, 
     *   coordMp, direct, T, charLength, props, 
     *   stateOld, 
c Write only -
     *   stateNew, field )
c
c
      include 'vaba_param.inc'
c
      dimension props(nprops),
     *          jElemUid(nblock), coordMp(nblock, *), 
     *          direct(nblock, 3, 3), T(nblock,3,3), 
     *          charLength(nblock),
     *          stateOld(nblock, nstatev), 
     *          stateNew(nblock, nstatev),
     *          field(nblock, nfieldv)
      character*80 cmname
c
      parameter( nrData=6 )
c     Input fracture model parameters
      parameter( a = 1.d0, b = 1.d0, r = 1.d0)

      
      
c Properties array
      character*3 cData(maxblk*nrData),qData(maxblk)
      dimension jData(maxblk*nrData),rData(maxblk,nrData)
      dimension eps(maxblk),pData(maxblk)
      real epsOld,damageOld,deltaeps,dam,damageNew
c  
c
      jStatus = 1
      call vgetvrm( 'S', rData, jData, cData, jStatus ) 
      call setField( nblock, nstatev, nfieldv, nrData,
     1  rData, stateOld, stateNew,  field) 
c
     	jStatus = 1
      call vgetvrm( 'PEEQ', eps, pData, qData, jStatus )

      
      
      do k = 1, nblock
c     Read value from the last increment
c     state 1 = Lode angle parameter
c     state 2 = Stress triaxiality
c     state 3 = Void growth rate
c     state 4 = plastic strain   
c     state 5 = Damage index
c     state 6 = Element deletion, 1 for no failure, 0 for failure
          
c     Plastic strain
        epsOld = stateOld(k,4)
c     Damage index
        damageOld = stateOld(k,5)
c     Calculationo f 
        deltaeps = eps(k) - epsOld
c     Void growth rate
        vg=stateOld(k,3)
c     D_new=D_old+void growth rate*plastic strian increment
        damageNew = damageOld + deltaeps*vg
c Element Deletion
        if(damageNew .gt. a) then 

            stateNew(k,6) = 0.d0
        endif
        stateNew(k,4) = eps(k)
        stateNew(k,5) = damageNew
      end do  
c
      return
      end
c
      subroutine setField(nblock, nstatev, nfieldv, nrdata,
     1  stress, stateOld,stateNew,field)
c
      include 'vaba_param.inc'
c
      dimension stateOld(nblock,nstatev),
     1  stateNew(nblock,nstatev),
     1  field(nblock,nfieldv),stress(nblock,nrData)
      real s11,s22,s33,s12,s23,s31,smises,sh123,smise
      do k = 1, nblock
c     Extract stress tensor
      	s11=stress(k,1)
      	s22=stress(k,2)
      	s33=stress(k,3)
      	s12=stress(k,4)
      	s23=stress(k,5)
      	s31=stress(k,6)
      	smises1=(s11-s22)**2.d0+(s22-s33)**2.d0+(s33-s11)**2.d0
      	smises2=(s12**2.d0+s23**2.d0+s31**2.d0)
      	smises3=0.5d0*(smises1+6.0d0*smises2)
      	smises4=sqrt(smises3)
          smises =smises4+1.d0
      	sh123=(s11+s22+s33)/3.0d0
          sig1=s11-sh123
          sig2=s22-sh123
          sig3=s33-sh123
c     stress triaxiality
          tri = sh123/(smises + 1.d0)
          J31=sig1*sig2*sig3+2*s12*s23*s31
          J32=sig1*(s23**2.d0)+sig2*(s31**2.d0)+sig3*(s12**2.d0)
          J3=J31-J32
          xi1=J3/(smises**3.d0+1.d0)
c     Lode angle parameter
          xi=27.0d0*xi1/2.d0
c     Effect of trixiality         
          dam1=2.718d0**(b*tri)
c     Effect of Lode angle parameter     
          dam2=2.718d0**(-r*xi)
          dam=dam1*dam2
        stateNew(k,1) = xi
        stateNew(k,2) = tri
        stateNew(k,3) = dam
      end do      
c
      return
      end