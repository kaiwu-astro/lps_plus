*     include file for fLPS.f

      INTEGER ::
     &           NHostStar,
     &           FindIndexByName(NMAX),
     &           HostStarNameList(NMAX),
     &           NPerturber, NActualPtb, LastNActualPtb(NMAX),
     &           LatestNBList(LMAX,NMAX),
     &           CurrentNBList(LMAX),
     &           LatestNBNumber,
     &           CurrentNBNumber,
     &           LatestPerturbersList(LMAX,NMAX),
     &           CurrentPerturbersList(LMAX),
     &           IndexStar, StarDevice, DiffList(LMAX),
     &           EscapedList(NMAX), Nescaped,
     &           tmpcount1, tmpcount2

      REAL*8 ::  LatestHostX1(NMAX),
     &           CurrentHostX1,
     &           LatestPerturbersX1(LMAX,NMAX),
     &           CurrentPerturbersX1(LMAX),
     &           RedoSearchTimeSep, NextRedoSearchTime,
     &           HostStarMassPref

      LOGICAL*1  IsNBListChanged, IsNBListUpdated, DoOutput,
     &           IsPerturberChanged, IsHostStarEscapedList(NMAX),
     &           IsHostStarXChanged, IsPerturberXChanged,
     &           IsPerturberOrXChanged


      COMMON/LPS/ NHostStar, FindIndexByName, HostStarNameList,
     & NPerturber, NActualPtb,
     & LatestNBList, CurrentNBList, LatestNBNumber, CurrentNBNumber,
     & LatestPerturbersList, CurrentPerturbersList, IndexStar,
     & StarDevice, DiffList, tmpcount1, tmpcount2, LastNActualPtb,
     & EscapedList, Nescaped,
     &
     & LatestHostX1,  CurrentHostX1,
     & LatestPerturbersX1, CurrentPerturbersX1, RedoSearchTimeSep,
     & NextRedoSearchTime, HostStarMassPref,
     &
     & IsNBListChanged,
     & IsNBListUpdated, DoOutput, IsPerturberChanged,
     & IsHostStarEscapedList, IsHostStarXChanged, IsPerturberXChanged,
     & IsPerturberOrXChanged
*----------------------------------------------------------------------