* fLPS.f : fortran implementation of LonelyPlanets Scheme
*     Kai Wu (XJTLU) created at 2020 Sept 23
*
*     LonelyPlanets is usually referred to as LPS, and sometimes LLP
*
*     Meaning of some variables used here: (Others can be easily infered with variable name)
*         HostStarNameList: names of host star. Names are the unique identifier of stars in Nbody6++GPU. Unsorted list. 
*         FindIndexByName: index of host star. Index are used among all lists in NBODY6++GPU, to reference mass, position, velocity, etc. Usage: IndexStar = FindIndexByName(HostStarNameList(i))
*         IsHostStarEscapedList: logical array. value = .true. or .false.
*         NHostStar: number of host stars
*         HostStarMassPref: prefered host star mass [solar mass]
*
*
*     Headers for host stars: 
*         Time[NB] Index  Name  Mass[Mdot, Solar mass]  
*         X1[pc]  X2  X3  V1[km/s]  V2  V3
*         A1[NBODY]  A2  A3  
*         J1[NBODY](Jerk, time derivative of acceleration) J2 J3
*         RS[Rdot, solar radius](Star radius)
*         KSTAR (Stellar type)

*     StarDevice: index of star in host star array, also deviceID of output file
*     IndexStar: index of star in Nbody6++GPU's BODY, which Nbody6 uses to identify stars. NOT NAME (which is NAME(IndexStar))
*     StarDevice一是用来索引HostStarIndex和IsHostStarEscapedList
*     二是用作fortran输出的设备号（本来想直接用host的name作为设备号，但设备号不能太大，而到后面1 Million particle的时候就不行了）
*     所以只会出现在这两个数组的下标和文件操作语句的设备号里，其他地方都应该用IndexStar

*     为了避免和主程序的循环变量I, J等有冲突，本程序所有循环都使用局域的Ind作为循环下标

*     因为有时候NBList长度很短，只有4，所以导致输出的ptblist有很远的，但并不影响任何计算

* 20230708: 发现不应该用ind而该用name追踪star，打补丁。建立简单的hash表来通过name查ind，每次输出时更新hash表。同时补上了之前缺失的判定escaper部分，现在经基本测试可以从escape.F里实时获取逃逸事件，并立即更新IsHostStarEscapedList（同时输出到LPSdiag.txt）
*---------------------------------------------------------------------
  
      SUBROUTINE LPSOutput
  
      include 'common6.h'
      INTEGER :: NameStar
  
      call Update_FindIndexByNameList
  
      DO StarDevice = 1, NHostStar 
         NameStar = HostStarNameList(StarDevice)
         IF (IsHostStarEscapedList(NameStar).EQV..true.) CYCLE ! if already escaped, skip
         IndexStar = FindIndexByName(NameStar)
         DoOutput = .false.
         CurrentNBNumber = LIST(1,IndexStar)
         NActualPtb = MIN(CurrentNBNumber,NPerturber)
  
         IF ((TTOT.GE.NextRedoSearchTime) .OR. 
     &       (LastNActualPtb(StarDevice).NE.NActualPtb)) THEN
            
            LastNActualPtb(StarDevice) = NActualPtb
            
            ! Step1: Get new NBList in the sequence of distance
            CurrentNBList(1:CurrentNBNumber) = 
     &               LIST(2:CurrentNBNumber+1,IndexStar)
            CALL quicksort_X(CurrentNBList,1,CurrentNBNumber,IndexStar)
  
            ! Step2: Check if perterbur list changes
            CurrentPerturbersList(1:NPerturber) = 
     &              CurrentNBList(1:NPerturber)
            CALL Compare_one_by_one_PList
            IF (IsPerturberChanged) THEN
               CALL Update_PerterberList
               DoOutput = .true.
            END IF
         END IF
  
         ! Step3: Check if host star position and ptb position changes
         IF (DoOutput.neqv..true.) THEN
            Call Host_star_position_changes ! update as
            IF (IsHostStarXChanged) THEN
               DoOutput = .true.
            ELSE
               Call Perturber_position_changes
               IF (IsPerturberXChanged) THEN
                  CALL Update_PerterberXList
                  DoOutput = .true.
               END IF
            END IF
         END IF
  
         IF (DoOutput) CALL LPS_Write_to_file
  
      END DO
  
      IF (TTOT.GE.NextRedoSearchTime) THEN ! Move forward the checkpoint
         NextRedoSearchTime = NextRedoSearchTime + RedoSearchTimeSep 
      END IF
  
      RETURN
      END
  
  
*---------------------------------------------------------------------
      SUBROUTINE LPSinit
  
*     Every host star uses 1 file
*     All neighbors of a host star use 1 file
*     Since LPS may use lots of units, and in order not to cause 
*       inconvenience to following code editors, 
*       I take large unit numbers starting from 10000 for LPS
  
      include 'common6.h'
  
      character( len = 7 ) :: IndChar
      INTEGER :: Ind
      
      DATA IsHostStarEscapedList/NMAX*.FALSE./
      tmpcount1 = 0
      tmpcount2 = 0
  
      NHostStar=100
      HostStarMassPref = 1.0
      CALL MakeHostStarList
      CALL Update_FindIndexByNameList
  
      NPerturber=10
      LastNActualPtb(1:NMAX) = NPerturber
      LatestNBList(1:LMAX,1:NMAX)=1
      CurrentNBList(1:LMAX)=2
      LatestNBNumber=0
      CurrentNBNumber=0
      LatestPerturbersList(1:LMAX,1:NMAX)=3
      CurrentPerturbersList(1:LMAX)=4
      
      LatestHostX1(NMAX)=20200925.15
      LatestPerturbersX1(LMAX,NMAX)=20201029.11
      
      RedoSearchTimeSep = 0.001 !in NB time ! 后记：如果设为0.001，则成了一个最大值，实际输出间隔基本完全由neighbor list变动决定，结果中位数一般是100yr以内
      NextRedoSearchTime = 0.0
  
      Nescaped=0
  
      DO Ind = 1, NHostStar
         WRITE(IndChar,'(i7)') HostStarNameList(Ind) !convert int to char
         OPEN(UNIT=10000+Ind,
     &        FILE='LPS_Host_'//trim(adjustl(IndChar))//'.txt',
     &        STATUS='UNKNOWN',FORM='FORMATTED',ACCESS='APPEND')
         OPEN(UNIT=20000+Ind,
     &        FILE='LPS_Perturber_'//trim(adjustl(IndChar))//'.txt',
     &        STATUS='UNKNOWN',FORM='FORMATTED',ACCESS='APPEND')
      END DO
      ! 如果后面要修改输出文件名，请注意不要在文件名中加入dot(.)。也就是，不要是用LPS.Host.123.txt这样的文件名。这是因为后面数据分析时，为了压缩输出文件，会采用.pkl.zst的后缀，而wukai现有的数据分析脚本，会把文件名用python:os.path.basename().split('.')[0]来获取，因此多于的点会导致文件名被错误识别。
  
      OPEN(UNIT=10000,
     &     FILE='LPSdiag.txt',
     &     STATUS='UNKNOWN',FORM='FORMATTED',ACCESS='APPEND')
  
      WRITE(10000,*) 
     &    "host_star_numbers", NHostStar, NEW_LINE("A"),
     &    "perturber_numbers", NPerturber, NEW_LINE("A"),
     &    "Mscale(to_Mdot)", ZMBAR, NEW_LINE("A"),
     &    "Rscale(to_pc)", RBAR, NEW_LINE("A"),
     &    "Vscale(to_kmps)", VSTAR, NEW_LINE("A"),
     &    "Tscale(to_Myr)", TSTAR
      call flush(10000) 
      ! WRITE (10000,*) "init: ",CurrentPerturbersList(1:NPerturber)
      
      END
*---------------------------------------------------------------------
      subroutine RecordEscaper(NAMEI)
      include 'common6.h'
      INTEGER :: NAMEI
      Nescaped = Nescaped + 1
      EscapedList(Nescaped) = NAMEI ! fortran数组下标默认从1开始，所以必须先+1再记录
      call Update_HostStarEscapeList
      END
! *---------------------------------------------------------------------
!       subroutine CheckHostStarEscaped(StarDeviceTMP)
!       include 'common6.h'
!       INTEGER :: Ind, StarDeviceTMP
!       ! check if host star's name is in EscapedList
!         ! if yes, set IsHostStarEscapedList(StarDevice) = .true.
!       DO Ind = 1, Nescaped
!          IF (HostStarNameList(StarDeviceTMP).EQ.EscapedList(Ind)) THEN
!             IsHostStarEscapedList(StarDeviceTMP) = .true.
!             ! write to LPSdiag that this host star escaped
!             WRITE(10000,*) "[ESCAPE] Host star ", HostStarNameList(StarDeviceTMP), 
!      &                     " escaped at ", TTOT*TSCALE, " Myr"
!             RETURN
!          END IF
!       END DO
!       END
*---------------------------------------------------------------------
      subroutine Update_HostStarEscapeList
      include 'common6.h'
      INTEGER :: Ind, NAME_NEW_ESCAPER
      NAME_NEW_ESCAPER = EscapedList(Nescaped)
      ! if this new escaping star is in HostStarNameList, write to LPSdiag
      DO Ind = 1, NHostStar
        IF (HostStarNameList(Ind).EQ.NAME_NEW_ESCAPER) THEN
            IsHostStarEscapedList(NAME_NEW_ESCAPER) = .true. 
            WRITE(10000,*) "[ESCAPE] Host star ", 
     &                     HostStarNameList(Ind), 
     &                     " escaped at ", TTOT*TSCALE, " Myr"
            call flush(10000)
        END IF
      END DO
      END
*---------------------------------------------------------------------
      subroutine Update_FindIndexByNameList
      include 'common6.h'
      INTEGER :: Ind, NameStartmp
      DO Ind = 1, NTOT
         NameStartmp = NAME(Ind)
         IF (NameStartmp.GT.0) FindIndexByName(NameStartmp) = Ind
      END DO
      END
*---------------------------------------------------------------------
      subroutine MakeHostStarList
      ! choose host star whose mass is the closest to the given mass
      include 'common6.h'
      INTEGER :: Ind, IndexCopy(NMAX), star_index ! make a copy, to be modified by sort
      REAL*8 :: MassDistanceToPref(NMAX), pref_in_nbmass
  
      pref_in_nbmass = HostStarMassPref / ZMBAR
      
      MassDistanceToPref(1:NTOT) = ABS(BODY(1:NTOT) - pref_in_nbmass)
  
      DO Ind = 1, NTOT
         IndexCopy(Ind) = Ind
      End DO
  
      call quicksort(MassDistanceToPref, 1, NTOT, IndexCopy)
      
      DO Ind = 1, NHostStar
         star_index = IndexCopy(Ind)
         HostStarNameList(Ind) = NAME(star_index)
      End DO
  
      END
*---------------------------------------------------------------------
      subroutine Update_PerterberList
      include 'common6.h'
  
!       if (NActualPtb.LT.NPerturber) THEN
!          LatestPerturbersList(1:NPerturber,IndexStar) = -2
!       end if 
! *如果NActualPtb < NPtb，就是这个neighbor list的长度太短
! *为了避免问题，先清空一遍PtbList，赋值-2方便debug
! *正常情况下，-2不应该出现在输出里，因为在输出函数LPS_Write里，空位是用-1填充的

! *上述改动弃用，已回滚，因为远的ptb不影响，这里改了后面全部要重写，太麻烦
!       LatestPerturbersList(1:NActualPtb,IndexStar) =
!      &         CurrentPerturbersList(1:NActualPtb)
      LatestPerturbersList(1:NPerturber,IndexStar) =
     &         CurrentPerturbersList(1:NPerturber)
  
      END
  
*---------------------------------------------------------------------
      subroutine Compare_one_by_one_PList
      include 'common6.h'
      INTEGER :: Ind
  
      IsPerturberChanged = .FALSE.
      DO Ind = 1, NPerturber
         IF (CurrentPerturbersList(Ind)
     &       .NE.LatestPerturbersList(Ind, IndexStar)) THEN
            IsPerturberChanged = .TRUE.
            EXIT
         END IF
      END DO
    !   IF (IsPerturberChanged) THEN
    !      tmpcount1 = tmpcount1 + 1
    !   ELSE
    !      tmpcount2 = tmpcount2 + 1
    !   END IF 
    !   WRITE(10000,*) tmpcount1, tmpcount2
      END
*---------------------------------------------------------------------
      subroutine Perturber_position_changes
      include "common6.h"
      INTEGER :: PerturberInd
      INTEGER :: Ind
         
      ! WRITE (10000,*) CurrentPerturbersList(1:NActualPtb)
      IsPerturberXChanged = .false.
      DO Ind = 1, NPerturber
         PerturberInd = CurrentPerturbersList(Ind)
         IF (X(1,PerturberInd).NE.
     &       LatestPerturbersX1(Ind,IndexStar)) THEN
            IsPerturberXChanged = .true.
            RETURN
         END IF
      END DO
  
      END
*---------------------------------------------------------------------
      subroutine Update_PerterberXList
      include "common6.h"
      INTEGER :: PerturberInd
      INTEGER :: Ind
      
      ! WRITE (10000,*) CurrentPerturbersList(1:NPerturber)
      DO Ind = 1, NPerturber
         PerturberInd = CurrentPerturbersList(Ind)
         LatestPerturbersX1(Ind,IndexStar) = X(1,PerturberInd)
      END DO
  
      END
*---------------------------------------------------------------------
      subroutine Host_star_position_changes
      include "common6.h"
  
      IsHostStarXChanged = .false.
      IF (X(1,IndexStar).NE.LatestHostX1(IndexStar)) THEN
         LatestHostX1(IndexStar) = X(1,IndexStar)
         IsHostStarXChanged = .true.
      END IF
  
      END
*---------------------------------------------------------------------
      subroutine LPS_Write_to_file
      include "common6.h"
      INTEGER :: PerturberInd
      INTEGER :: Ind
  
      !Write host star data
      WRITE (10000+StarDevice,*) 
     &   TTOT*TSCALE, IndexStar, NAME(IndexStar),
     &   BODY(IndexStar)*ZMBAR,
     &   X(1:3,IndexStar)*RBAR, XDOT(1:3,IndexStar)*VSTAR
*     &   F(1:3,IndexStar), 
*     &   FDOT(1:3,IndexStar),
*     &   RS(IndexStar)*RBAR, ! RS: neighbor radius [NB], not star radius in manual
*     &   KSTAR(IndexStar)
  
      !Write Perturber data
      DO Ind = 1, NPerturber
         PerturberInd = LatestPerturbersList(Ind, IndexStar)
         WRITE (20000+StarDevice,*)
     &      TTOT*TSCALE, PerturberInd, NAME(PerturberInd),
     &      BODY(PerturberInd)*ZMBAR,
     &      X(1:3,PerturberInd)*RBAR, XDOT(1:3,PerturberInd)*VSTAR
*     &      F(1:3,PerturberInd), 
*     &      FDOT(1:3,PerturberInd),
*     &      RS(PerturberInd)*RBAR,
*     &      KSTAR(PerturberInd)
      END DO
      
*    Neighbor List长度小于需要的Perturber数量时输出诊断信息
    !   IF (NActualPtb.LT.NPerturber) THEN
    !      WRITE(10000,*) "NNB=",NActualPtb," < NPtb=",NPerturber,
    !  &    " at host", IndexStar
    !   END IF
      
      END
*---------------------------------------------------------------------
      subroutine QuickSort_X(arr, low, high, centerind)
      include "common6.h"
      REAL*8 distsq(LMAX)
      INTEGER particleind, ind, arr(*), low, high, centerind
      do ind = low, high
         particleind = arr(ind)
         distsq(ind) = (X(1,particleind) - X(1,centerind))**2 + 
     &                 (X(2,particleind) - X(2,centerind))**2 +
     &                 (X(3,particleind) - X(3,centerind))**2
      end do
      call quicksort(distsq,low,high,arr)
  
      END
*---------------------------------------------------------------------
      subroutine quicksort(array,lowin,highin,witharr)
  
      ! array: Sort by
      ! lowin & highin: range to sort of array
      ! witharr: the array whose element sequence changes 
      !          with the sort-by array, which is usually the index
      
      ! This is the non-recursive quicksort subroutine
      ! originally from https://www.mjr19.org.uk/IT/sorts/sorts.f90
      ! modified by Kai Wu for sorting with moving index
  
      real*8 :: array(*)
      real*8 :: temp,pivot
      integer :: i,j,left,right,low,high,lowin,highin
      integer :: witharr(*), withtmp
      integer :: stack(2,100),stack_ptr
      
      ! low=1
      ! high=size(array)
      low = lowin
      high = highin
      stack_ptr=1
      do
         if (high-low.lt.50) then ! use insertion sort on small arrays
            do i=low+1,high
               temp=array(i)
               withtmp = witharr(i)
               do j=i-1,low,-1
                  if (array(j).le.temp) exit
                  array(j+1)=array(j)
                  witharr(j+1) = witharr(j)
               enddo
               array(j+1)=temp
               witharr(j+1) = withtmp
            enddo
            ! now pop from stack
            if (stack_ptr.eq.1) return
            stack_ptr=stack_ptr-1
            low=stack(1,stack_ptr)
            high=stack(2,stack_ptr)
            cycle
         endif
         
         ! find median of three pivot
         ! and place sentinels at first and last elements
         temp=array((low+high)/2)
         array((low+high)/2)=array(low+1)
         withtmp = witharr((low+high)/2)
         witharr((low+high)/2) = witharr(low+1)
         if (temp.gt.array(high)) then
            array(low+1)=array(high)
            array(high)=temp
            witharr(low+1)=witharr(high)
            witharr(high)=withtmp
         else
            array(low+1)=temp
            witharr(low+1)=withtmp
         endif
         if (array(low).gt.array(high)) then
            temp=array(low)
            array(low)=array(high)
            array(high)=temp
            withtmp=witharr(low)
            witharr(low)=witharr(high)
            witharr(high)=withtmp
         endif
         if (array(low).gt.array(low+1)) then
            temp=array(low)
            array(low)=array(low+1)
            array(low+1)=temp
            withtmp=witharr(low)
            witharr(low)=witharr(low+1)
            witharr(low+1)=withtmp
         endif
         pivot=array(low+1)
         
         left=low+2
         right=high-1
         do
            do while(array(left).lt.pivot)
               left=left+1
            enddo
            do while(array(right).gt.pivot)
               right=right-1
            enddo
            if (left.ge.right) exit
            temp=array(left)
            array(left)=array(right)
            array(right)=temp
            withtmp=witharr(left)
            witharr(left)=witharr(right)
            witharr(right)=withtmp
            left=left+1
            right=right-1
         enddo
         if (left.eq.right) left=left+1
         if (left.lt.(low+high)/2) then
            stack(1,stack_ptr)=left
            stack(2,stack_ptr)=high
            stack_ptr=stack_ptr+1
            high=left-1
         else
            stack(1,stack_ptr)=low
            stack(2,stack_ptr)=left-1
            stack_ptr=stack_ptr+1
            low=left
         endif
      enddo
      end 
  
