      program help
c**************************************************************************
c
c     designed to create inputs and aid in keeping track of molecule
c     groups during the creation of gcmc simulations, dlpoly simulations
c     and cpmd simulations of a framework.
c
c**************************************************************************

      implicit none
      integer, parameter :: lenrec=80
      integer, parameter :: nturd=205
      integer, parameter :: npdb=206

      real(8), parameter :: angs2bohr=1.889725989d0
      logical init,lsym,lcpmd,geom,elpot,cpmdread,lcpmdord
      logical reptwrite,reptread,fldwrite,cnfgwrite,lvdw,lpdb

      integer natms,nmols,ntpsit,i,j,ntpatm,ind,ls,ntpmol,na,nb,nc,nvdw
      integer, allocatable :: nsym(:),molind(:),ltype(:),printorder(:)
      integer, allocatable :: cpmdorder(:),natm(:),latm(:,:)
      integer, allocatable :: lstsym(:,:),symorder(:),lstmol(:)
      integer, allocatable :: numatoms(:)
      real(8), allocatable :: xxx(:), yyy(:), zzz(:)
      real(8), allocatable :: charges(:),mass(:),eps(:),sig(:)
      character*1 title(lenrec)
      character*1, allocatable :: rootname(:)
      character*2, allocatable :: atom(:),unqatm(:)
      character*3, allocatable :: atmname(:),unqsit(:),vdwatm(:)
      real(8), dimension(9) :: cell,rcell
      real(8), dimension(10) :: celprp
      real(8) alpha,beta,gamma,alen,blen,clen,chg

      character*1 record(lenrec)
      character*70 pdbfile,ljfile,cpmdfile

      save record,xxx,yyy,zzz,charges,mass,atom,atmname
      save unqsit,nsym,lstsym,molind,ltype,printorder,lstmol
      save cpmdorder,unqatm,latm,natm,rootname,symorder,numatoms
      save eps,sig,vdwatm
      data init,lcpmd,geom,elpot/.false.,.false.,.false.,.false./
      data cpmdread,lcpmdord,reptwrite/.false.,.false.,.false./
      data reptread,fldwrite,cnfgwrite/.false.,.false.,.false./
      data lvdw,lpdb/.false.,.false./
c     initial default case
      lsym=.true.
      chg=0.d0
      call argparse(pdbfile,init,lcpmd,cpmdread,reptwrite,
     &reptread,geom,elpot,lvdw,fldwrite,na,nb,nc,cnfgwrite,lpdb,
     &ljfile,lsym)

      if(init)then
c     need to look for existing TURDS file before opening this?
c         open(nturd,file='TURDS',status='new',form='unformatted')
         open(nturd,file='TURDS',form='unformatted')
         open(npdb,file=pdbfile,status='old')
c     scan pdb file to determine number of atoms, molecules to allocate
c     arrays
         call getroot(pdbfile,ind)
         call initscan(npdb,nmols,natms,ntpsit)
         allocate(xxx(natms),yyy(natms),zzz(natms),printorder(natms))
         allocate(atom(natms),atmname(natms),mass(natms))
         allocate(charges(natms),molind(natms),ltype(natms))
         allocate(lstsym(natms,natms),nsym(natms),unqsit(natms))
         allocate(latm(natms,natms),natm(natms),unqatm(natms))
         allocate(cpmdorder(natms),symorder(natms),numatoms(natms))
         allocate(lstmol(natms))
c     ** note, assumes pdb file is a representation of the unit cell **
         do i=1,natms
           charges(i)=0.d0
           ltype(i)=i
           nsym(i)=1
           natm(i)=1
           cpmdorder(i)=i
           printorder(i)=i
           mass(i)=0.d0
           molind(i)=1
           ltype(i)=i
           symorder(i)=i
         enddo
         call readpdb(npdb,pdbfile,nmols,natms,lsym,ntpsit,
     & ntpatm,ntpmol,cell,alen,blen,clen,alpha,beta,gamma,title) 

      
         call getcell(cell,alen,blen,clen,alpha,beta,gamma) 



c     write all information to the TURDS file
         call writeturd
     &(nturd,natms,nmols,ntpsit,ntpatm,ntpmol,alen,blen,clen,
     &alpha,beta,gamma,cell,lsym,lcpmdord,ind,lvdw,nvdw,title)

         close(npdb)
      endif

      if(lvdw)then
        lvdw=.false.
        if(.not.allocated(xxx))call readturd
     &(nturd,natms,nmols,ntpsit,ntpatm,ntpmol,alen,blen,clen,
     &alpha,beta,gamma,cell,lsym,lcpmdord,ind,lvdw,nvdw,title)
        lvdw=.true.
        allocate(vdwatm(natms),eps(natms),sig(natms))
        call vdwread(nvdw,ljfile)

c     write all information to the TURDS file
        call writeturd
     &(nturd,natms,nmols,ntpsit,ntpatm,ntpmol,alen,blen,clen,
     &alpha,beta,gamma,cell,lsym,lcpmdord,ind,lvdw,nvdw,title)

      endif 

      if(lpdb)then
        if(.not.allocated(xxx))call readturd
     &(nturd,natms,nmols,ntpsit,ntpatm,ntpmol,alen,blen,clen,
     &alpha,beta,gamma,cell,lsym,lcpmdord,ind,lvdw,nvdw,title)
c     adjust the xyz coordinates to have the center of the cell as the 
c     origin.  molecules whose centroids are within the cell will remain
c     others will be shifted by the vectors.
         
c        call images(natms,cell,xxx,yyy,zzz)
        call writepdb(ind,title,natms,alen,blen,clen,alpha,beta,
     &gamma)

      endif
      if(cpmdread)then
        
        if(.not.allocated(xxx))call readturd
     &(nturd,natms,nmols,ntpsit,ntpatm,ntpmol,alen,blen,clen,
     &alpha,beta,gamma,cell,lsym,lcpmdord,ind,lvdw,nvdw,title)

        call readcpmd(ind,geom,elpot,natms,ntpatm)

c     write all information to the TURDS file
        call writeturd
     &(nturd,natms,nmols,ntpsit,ntpatm,ntpmol,alen,blen,clen,
     &alpha,beta,gamma,cell,lsym,lcpmdord,ind,lvdw,nvdw,title)
      endif

      if(cnfgwrite)then
        
        if(.not.allocated(xxx))call readturd
     &(nturd,natms,nmols,ntpsit,ntpatm,ntpmol,alen,blen,clen,
     &alpha,beta,gamma,cell,lsym,lcpmdord,ind,lvdw,nvdw,title)
      
        call writeconfig(natms,ntpmol,ntpsit,title,
     &na,nb,nc,cell)
      endif 

      if(fldwrite)then
        
        if(.not.allocated(xxx))call readturd
     &(nturd,natms,nmols,ntpsit,ntpatm,ntpmol,alen,blen,clen,
     &alpha,beta,gamma,cell,lsym,lcpmdord,ind,lvdw,nvdw,title)
      
        call writefield(natms,ntpmol,ntpsit,title,
     &na,nb,nc,lvdw,nvdw,lsym)
      endif 
      if(reptread)then
        if(.not.allocated(xxx))call readturd
     &(nturd,natms,nmols,ntpsit,ntpatm,ntpmol,alen,blen,clen,
     &alpha,beta,gamma,cell,lsym,lcpmdord,ind,lvdw,nvdw,title)
      
        if(.not.lcpmdord)then
           write(*,"('ERROR - no CPMD info found in TURDS file. 
     &Please make sure you write a cpmd input file before reading 
     &this file')")
           stop
        endif
        call readrepeat(ind,ntpsit,lsym)
c     write all information to the TURDS file

        call writeturd
     &(nturd,natms,nmols,ntpsit,ntpatm,ntpmol,alen,blen,clen,
     &alpha,beta,gamma,cell,lsym,lcpmdord,ind,lvdw,nvdw,title)
      endif


      if(reptwrite)then
        if(.not.allocated(xxx))call readturd
     &(nturd,natms,nmols,ntpsit,ntpatm,ntpmol,alen,blen,clen,
     &alpha,beta,gamma,cell,lsym,lcpmdord,ind,lvdw,nvdw,title)

       
        if(.not.lcpmdord)then
           write(*,"('ERROR - no CPMD info found in TURDS file. 
     &Please make sure you write a cpmd input file before writing 
     &this file')")
           stop
        endif
        call writerepeat(ind,ntpsit,lsym,chg)
      endif

      if(lcpmd)then
        
c       check to see if the arrays have already been allocated:
c       this would mean that the pdb file has just been read and
c       we do not need to populate arrays from TURDS
        if(.not.allocated(xxx))call readturd
     &(nturd,natms,nmols,ntpsit,ntpatm,ntpmol,alen,blen,clen,
     &alpha,beta,gamma,cell,lsym,lcpmdord,ind,lvdw,nvdw,title)
        call writecpmd(ind,lcpmdord,geom,elpot,natms,ntpatm,
     &cell,chg,lsym,title)
c     write all information to the TURDS file
        call writeturd
     &(nturd,natms,nmols,ntpsit,ntpatm,ntpmol,alen,blen,clen,
     &alpha,beta,gamma,cell,lsym,lcpmdord,ind,lvdw,nvdw,title)

      endif
      if(.not.allocated(xxx))call readturd
     &(nturd,natms,nmols,ntpsit,ntpatm,ntpmol,alen,blen,clen,
     &alpha,beta,gamma,cell,lsym,lcpmdord,ind,lvdw,nvdw,title)
      close(nturd)
      contains

      subroutine writepdb(ind,title,natms,alen,blen,clen,alpha,beta,
     &gamma)
c*********************************************************************
c     write a pdb file with the new order of the atoms 
c*********************************************************************
      implicit none
      integer natms,i,ind
      real(8), parameter :: rad2deg=57.2957795d0
      real(8) alen,blen,clen,alpha,beta,gamma
      character*1 title(lenrec)
      character(len=ind) root
      character(len=ind+8) pdbfile


      do i=1,ind
        root(i:i)=rootname(i)
      enddo

      pdbfile=root//"_new.pdb"

      open(77,file=pdbfile) 
      write(77,"('REMARK',3x,71a)")title(1:71)
      write(77,"('CRYST1',3x,f6.3,3x,f6.3,3x,f6.3,1x
     &,f6.2,1x,f6.2,1x,f6.2,1x,'P1')")alen,blen,clen,
     &alpha*rad2deg,beta*rad2deg,gamma*rad2deg
      do i=1,natms
        write(77,"('ATOM  ',i5,1x,a3,2x,'MOL ',i5,4x,3f8.3,
     &2x,f4.2,2x,f4.2,10x,a2)")i,atmname(i),molind(i),xxx(i),yyy(i),
     & zzz(i),1.d0,0.d0,atom(i)
      enddo
    
      write(77,"('TER')") 

      close(77)

      return
      end subroutine writepdb
      subroutine vdwread(nvdw,ljfile)
c*********************************************************************
c     read in lj paramter input file and mix the params
c*********************************************************************
      implicit none
      logical safe
      integer i,j,nvdw,idum
      character*8 atn
      character*70 ljfile

      safe=.true.
      open(77,file=ljfile)
c     initial parameter count
      nvdw=0
      do while(safe)
        call getrec(safe,77)
        call strip(record,lenrec)
        if(record(1).ne.' '.and.safe)then
          nvdw=nvdw+1
        endif
      enddo
      close(77)
      
      open(77,file=ljfile)
      do i=1,nvdw
        call getrec(safe,77)
        call strip(record,lenrec)
        call getword(atn,record,3,lenrec)
        vdwatm(i)=atn(1:3)
        eps(i)=dblstr(record,lenrec,idum)
        sig(i)=dblstr(record,lenrec,idum)
      enddo
      close(77)
      return
      end subroutine vdwread

      subroutine writeconfig(natms,ntpmol,ntpsit,title,
     &na,nb,nc,cell)
c*********************************************************************
c     writes a CONFIG file for dl_poly
c*********************************************************************
      implicit none
      integer natms,ntpmol,ntpsit,na,nb,nc,i,j,ia,ib,ic
      integer atmcount,cmol
      integer, dimension(natms) :: buffer
      character*1 title(lenrec)
      real(8), dimension(9) :: cell,ncell,tcell
      real(8) x,y,z

      ncell(1:3)=cell(1:3)*na
      ncell(4:6)=cell(4:6)*nb
      ncell(7:9)=cell(7:9)*nc

      open(66,file='CONFIG')
      write(66,'(80a)')title
      write(66,'(2i10)')0,3
      write(66,'(3f20.12)')ncell

      atmcount=0
      do i=1,ntpmol
        cmol=0
        do j=1,natms
          if(i.eq.molind(j))then
            atmcount=atmcount+1
c            write(66,'(a2,10x,i4,/,3f20.12)')
c     &atom(j),atmcount,xxx(j),yyy(j),zzz(j)
            cmol=cmol+1
            buffer(cmol)=j
          endif
        enddo
        do ia=0,na-1
          do ib=0,nb-1
            do ic=0,nc-1
              tcell(1:3)=cell(1:3)*ia
              tcell(4:6)=cell(4:6)*ib
              tcell(7:9)=cell(7:9)*ic
              do j=1,cmol
                atmcount=atmcount+1
                x=xxx(buffer(j))+tcell(1)+tcell(4)+tcell(7)
                y=yyy(buffer(j))+tcell(2)+tcell(5)+tcell(8)
                z=zzz(buffer(j))+tcell(3)+tcell(6)+tcell(9)
                write(66,'(a2,10x,i4,/,3f20.12)')
     &atom(buffer(j)),atmcount,x,y,z
                
              enddo
            enddo
          enddo
        enddo
      enddo
      

      return
      end subroutine writeconfig
      subroutine writefield(natms,ntpmol,ntpsit,title,
     &na,nb,nc,lvdw,nvdw,lsym)
c*********************************************************************
c     writes a field file for dl_poly
c*********************************************************************
      implicit none
      logical lvdw,lsym
      integer natms,ntpmol,ntpsit,i,j,na,nb,nc,atcount,nvdw,satm
      integer vcnt
      real(8) schg,ep,si
      character*1 title(lenrec)

      open(55,file='FIELD')
      write(55,'(80a)')title

      write(55,"('UNITS',15x,'kcal')")
      write(55,"('molecular types',i5)")ntpmol
      do i=1,ntpmol
        write(55,"('Molecule',i3)")i
        write(55,"('NUMMOLS',2x,i3)")na*nb*nc
        write(55,"('ATOMS',4x,i3)")numatoms(i)
        atcount=0
        do j=1,natms
          if(i.eq.molind(j))then
c     get the particular charge type for the atom
            if(lsym)then
              satm=ltype(j)
            else
              satm=j
            endif
            schg=charges(satm)
            write(55,'(a2,13x,f7.4,3x,f20.15,5x,i1,5x,i1)')atom(j),
     &mass(j),schg,1,1
            atcount=atcount+1
          endif
        enddo
        if(atcount.ne.numatoms(i))then
          write(*,"('ERROR - there should be ',i2,' atoms in molecule ',
     &i2,' but ',i2,' are reported here!')")numatoms(i),i,atcount
          stop
        else
          write(55,"('finish')")
        endif 
      enddo
      if(lvdw)then
        vcnt=0
        do i=0,nvdw
          vcnt=vcnt+nvdw-i
        enddo
        write(55,"('VDW',12x,i4)")vcnt
        do i=1,nvdw
          do j=i,nvdw
            if((eps(i).ne.0.d0).and.(eps(j).ne.0.d0))then
               ep=sqrt(eps(i)*eps(j))
            else
               ep=0.d0
            endif
            if((sig(i).ne.0.d0).and.(sig(j).ne.0.d0))then
               si=(sig(i)+sig(j))/2.d0
            else
               si=0.d0
            endif
            write(55,"(a3,7x,a3,5x,'lj',3x,f8.5,7x,f8.5)")
     &vdwatm(i),vdwatm(j),ep,si
          enddo
        enddo
      endif
      write(55,"('close')")
      close(55)
      return
      end subroutine writefield

      subroutine readrepeat(ind,ntpsit,lsym)
c*********************************************************************
c     read repeat input file ELPOT.esp_fit.out 
c     and assign charges to all the atoms 
c*********************************************************************
      implicit none
      logical lsym,done,safe
      integer i,j,jsite,ind,ntpsit,idum,natm
      real(8) ncharge
      character(len=ind) root
      character(len=ind+12) reptfile

      data done,safe/.false.,.true./

      do i=1,ind
        root(i:i)=rootname(i)
      enddo
      reptfile=root//'.esp_fit.out'

      open(10,file=reptfile,status='old')

      do while(.not.done)

        call getrec(safe,10)
        if(findstring('Fitted charges ordered',record,idum))
     &done=.true.
      enddo

      if(done)then
        do i=1,ntpsit
          call getrec(safe,10)
          do j=30,lenrec
            record(j-30)=record(j)
          enddo
          jsite=symorder(i) 
          natm=intstr(record,lenrec,idum)
           
          charges(jsite)=dblstr(record,lenrec,idum)
        enddo
      endif



      close(10)
      return
      end subroutine readrepeat
      subroutine writerepeat(ind,ntpsit,lsym,chg)
c*********************************************************************
c     write repeat input files, REPEAT_params.inp and 
c     connectivity.ff if lsym=.true.
c*********************************************************************
      implicit none

      integer i,j,ind,symmetry,ntpsit,iatm
      logical lsym
      real(8) chg
      character(len=ind) root
      character(len=ind+5) cubefile

      do i=1,ind
        root(i:i)=rootname(i)
      enddo

      cubefile=root//'.cube'
     
      if(lsym)then
         symmetry=1
      else
         symmetry=0
      endif


      open(9,file='REPEAT_param.inp')

      write(9,"('ESP file name',/,a)")cubefile

      write(9,"('Fit molecular(0) or periodic(1:default) system?',/,i1)"
     &)1
      write(9,"('van der Waals scaling factor (default = 1.0)',/,f4.2)")
     &1.d0
      write(9,"('Apply RESP penalties?, no(0:default), yes(1)',/,i1)")0

      write(9,"('Read cutoff radius? no(0), yes(1:default)',/,i1)")1

      write(9,"('If flag above=1 provide R_cutoff next (in
     &Bohrs)',/,f8.5)")20.d0

      write(9,"('Apply symmetry restrain? no(0:default), yes(1)',/,i1)")
     &symmetry

      write(9,"('Use Goddard-type restraints? no(0:default),
     &yes(1)',/,i1)")0

      write(9,"('If flag above=1 then provide weight next',/,f7.5)")
     &0.d0
      
      write(9,"('Enter total charge of the system',/,f7.5)")chg

      if(lsym)then
         open(10,file='connectivity.ff')
         write(10,'(i3)',advance='no')ntpsit
         do i=1,ntpsit
           write(10,"(/,'#',a3,/,i3)")unqsit(i),nsym(i)
           do j=1,nsym(i)
             iatm=lstsym(i,j)
             write(10,"(i3,1x)",advance='no')printorder(iatm)
             
           enddo

         enddo
      endif

      return
      end subroutine writerepeat
      subroutine getroot(pdbfile,ind)
c*********************************************************************
c     gather the pdb filename before the .pdb
c     this will enable the writing of other output files
c*********************************************************************
      implicit none
      integer i,ind
      character*70 pdbfile

      do i=1,70
        if(pdbfile(i:i).eq.'.')then
           ind=i-1
        endif
      enddo
      allocate(rootname(ind))
      do i=1,ind
        rootname(i)=pdbfile(i:i)
      enddo
      return
      end subroutine getroot

      subroutine readcpmd(ind,geom,elpot,natms,ntpatm)
c***********************************************************************
c     
c    write cpmd input file for submission 
c     
c***********************************************************************
      implicit none

      integer, parameter :: outfile=3
      real(8), parameter :: bohr2angs=0.529177249d0
      logical safe,done,atmread
     
      character(len=ind+13)cpmdfile
      character(len=ind) root
      logical geom,elpot
      integer ind,i,j,jatm,natms,ntpatm,icpmd,idum
      real(8) xcoord,ycoord,zcoord 


      data done,safe,atmread/.false.,.true.,.false./
      safe=.true.
      
      do i=1,ind
        root(i:i)=rootname(i)
      enddo

      cpmdfile=root//'_cpmd_geo.out'
      open(outfile,file=cpmdfile,status='old')

c     WARNING THIS FORMAT MIGHT ONLY BE APPLICABLE TO
c     CPMD VERSION 3.13.2
      do while(.not.done)
        call getrec(safe,outfile)
        if(findstring('***************',record,idum))then
          call getrec(safe,outfile)
          call getrec(safe,outfile)
          if(findstring('FINAL RESULTS',record,idum))then
            atmread=.true.
            done=.true.
            do i=1,4
              call getrec(safe,outfile)
            enddo
 
          endif
        endif
      enddo
      if(atmread)then
        do i=1,natms
          jatm=cpmdorder(i)
c         delete first 7 entries in line 
          call getrec(safe,outfile)
          do j=7,lenrec
            record(j-7)=record(j)
          enddo
          xxx(jatm)=dblstr(record,lenrec,idum)*bohr2angs
          yyy(jatm)=dblstr(record,lenrec,idum)*bohr2angs
          zzz(jatm)=dblstr(record,lenrec,idum)*bohr2angs
          
        enddo
      endif

      close(outfile)
      return
      end subroutine readcpmd
      subroutine writecpmd(ind,lcpmdord,geom,elpot,natms,ntpatm,cell,
     &chg,lsym,title)
c***********************************************************************
c     
c    write cpmd input file for submission 
c     
c***********************************************************************
      implicit none

      integer, parameter :: outfile=3
      character*1 title(80)
      character*8 func
      character*24 pseudo
      character(len=ind+12)cpmdfile
      character(len=ind) root
      logical geom,elpot,lcpmdord,chk,lsym
      integer ind,i,j,k,jatm,natms,ntpatm,icpmd,hstart,hend,sym,tp
      real(8) cut,chg,cell(9)

      data tp/0/
      do i=1,ind
        root(i:i)=rootname(i)
      enddo
      if(geom)then
        cpmdfile=root//'_cpmd_geo.in'
      elseif(elpot)then
        cpmdfile=root//'_cpmd_elp.in'
      endif
      
      open(outfile,file=cpmdfile)

c     cpmd parameters
      cut=60.d0

      func='PBE'
      pseudo='_VDB_PBE.psp  FORMATTED'

      hstart=0
      hend=0
c     write CPMD info
      write(outfile,"('&INFO',/,3x,80a,/,'&END')")title

c     CPMD section
      write(outfile,"('&CPMD')")
      if(geom)write(outfile,"(3x,'OPTIMIZE GEOMETRY XYZ',/,
     &3x,'CONVERGENCE GEOMETRY',/,6x,d7.1)")
     & 5.d-4
      if(elpot)write(outfile,"(3x,'ELECTROSTATIC POTENTIAL',/,
     &3x,'OPTIMIZE WAVEFUNCTION')")
      write(outfile,"(3x,'CONVERGENCE ORBITALS',/,6x,d7.1)")
     & 1.d-5
      write(outfile,"(3x,'PCG MINIMIZE',/,'&END',/)")

c     SYSTEM section of the input file      
      write(outfile,"('&SYSTEM',/,3x,'CELL VECTOR')")

      write(outfile,'(3f12.5)')cell*angs2bohr

      write(outfile,"(3x,'CUTOFF',/,6x,f4.1,/,3x,'CHARGE',/,6x,f4.1)")
     &cut,chg
      write(outfile,"('&END')")

c     DFT section
      write(outfile,"('&DFT',/,3x,'FUNCTIONAL ',a9,/,'&END')")func


c     ATOMS section
      write(outfile,"('&ATOMS')")
      icpmd=0
      do i=1,ntpatm
        call atmline(outfile,unqatm(i),pseudo)
        if((unqatm(i).eq.' H').or.(unqatm(i).eq.'H '))then
          write(outfile,"(3x,'LMAX=S',/,i4)")natm(i)
          hstart=icpmd
          hend=icpmd+natm(i)
        else
          write(outfile,"(3x,'LMAX=P',/,i4)")natm(i)
        endif
      
        do j=1,natm(i)
          icpmd=icpmd+1
          jatm=latm(i,j)
          
          lcpmdord=.true.
          cpmdorder(icpmd)=jatm
          printorder(jatm)=icpmd
          write(outfile,10)xxx(jatm)*angs2bohr,yyy(jatm)*angs2bohr,
     &zzz(jatm)*angs2bohr

c     populate symorder array, this makes the repeat charges
c     easy to assign to the correct atom type
          if(lsym)then
            sym=ltype(jatm)

            chk=.true.
            do k=1,tp
              if(sym.eq.symorder(k))chk=.false.
            enddo
            if(chk)then
              tp=tp+1
              symorder(tp)=sym
            endif
          else
            symorder=cpmdorder
          endif


        enddo
      enddo 

      if(geom)then
        write(outfile,"(3x,'&CONSTRAINTS')")
        write(outfile,"(6x,'&FIX SEQUENCE',
     &/,9x,i3,2x,i3)")1,hstart
        if(hend.ne.natms)then
          write(outfile,"(6x,'&FIX SEQUENCE',/,
     &9x,i3,2x,i3)")hend+1,natms
        endif
        write(outfile,"(3x,'&END CONSTRAINTS')")
      endif
      write(outfile,"('&END')")

10    format(2x,3f12.6)
      close(outfile)
      return
      end subroutine writecpmd 


      subroutine atmline(outfile,atmname,pseudo)
c*********************************************************************
c     quick subroutine to solve the problem of elements containing
c     1 or 2 characters eg. Zn vs. H
c*********************************************************************
      implicit none
      integer outfile
      character*24 pseudo
      character*2 atmname

      if(atmname(1:1)//' '.eq.atmname)then
        write(outfile,"('*',a1,a24)")atmname(1:1),pseudo
      elseif(' '//atmname(2:2).eq.atmname)then
        write(outfile,"('*',a1,a24)")atmname(2:2),pseudo
      else
        write(outfile,"('*',a2,a24)")atmname,pseudo
      endif
      return
      end subroutine atmline
      subroutine writeturd
     &(nturd,natms,nmols,ntpsit,ntpatm,ntpmol,alen,blen,clen,
     &alpha,beta,gamma,cell,lsym,lcpmdord,ind,lvdw,nvdw,title)
c***********************************************************************
c     
c    write unformatted turd file, this will contain all information
c    about the framework necessary to write files etc..
c     
c***********************************************************************
      implicit none
      logical lsym,lcpmdord,lvdw
      character*1 title(80)
      integer nturd,natms,nmols,ntpsit,i,j,ntpatm,ind,ntpmol,nvdw
      real(8) alen,blen,clen,alpha,beta,gamma
      real(8), dimension(9) :: cell

      open(nturd,file='TURDS',status='unknown',form='unformatted')
c     write integer values first
      write(nturd)natms,nmols,ntpsit,ntpatm,ntpmol,ind

c     real values.. not sure i need these any more
      write(nturd)alen,blen,clen,alpha,beta,gamma

c     write cell dimensions
      write(nturd)(cell(i),i=1,9)

c     write atomistic properties
      do i=1,natms
        write(nturd)xxx(i)
        write(nturd)yyy(i)
        write(nturd)zzz(i)
        write(nturd)atom(i)
        write(nturd)atmname(i)
        write(nturd)mass(i)
        write(nturd)molind(i)
        write(nturd)ltype(i)
      enddo
      
      write(nturd)(numatoms(i),i=1,ntpmol)

c     write symmetry related stuff
      write(nturd)lsym
 
      do i=1,ntpsit
        write(nturd)unqsit(i)
        write(nturd)nsym(i)
        write(nturd)(lstsym(i,j),j=1,nsym(i))
      enddo

      do i=1,ntpatm
        write(nturd)unqatm(i)
        write(nturd)natm(i)
        write(nturd)(latm(i,j),j=1,natm(i))
      enddo

      write(nturd)(rootname(i),i=1,ind)
      write(nturd)(title(i),i=1,80)
      write(nturd)(symorder(i),i=1,ntpsit)

      write(nturd)(charges(i),i=1,ntpsit)

      write(nturd)lcpmdord
      if(lcpmdord)write(nturd)(cpmdorder(i),i=1,natms)
      if(lcpmdord)write(nturd)(printorder(i),i=1,natms)

      write(nturd)lvdw
      if(lvdw)then
        write(nturd)nvdw
        write(nturd)(vdwatm(i),i=1,nvdw)
        write(nturd)(eps(i),i=1,nvdw)
        write(nturd)(sig(i),i=1,nvdw)
      endif 


      return
      end subroutine writeturd


      subroutine readturd
     &(nturd,natms,nmols,ntpsit,ntpatm,ntpmol,alen,blen,clen,
     &alpha,beta,gamma,cell,lsym,lcpmdord,ind,lvdw,nvdw,title)
c***********************************************************************
c     
c    read unformatted turd file, this will contain all information
c    about the framework necessary to write files etc..
c     
c***********************************************************************
      implicit none
      character*1 title(80)
      logical lsym,lcpmdord,lvdw
      integer nturd,natms,nmols,ntpsit,ntpatm,i,j,ind,ntpmol,nvdw
      real(8) alen,blen,clen,alpha,beta,gamma
      real(8), dimension(9) :: cell


      open(nturd,file='TURDS',status='old',form='unformatted')
c     read integer values first
      read(nturd)natms,nmols,ntpsit,ntpatm,ntpmol,ind

c     allocate arrays for reading further info
      allocate(xxx(natms),yyy(natms),zzz(natms),printorder(natms))
      allocate(atom(natms),atmname(natms),mass(natms))
      allocate(charges(ntpsit),molind(natms),ltype(natms))
      allocate(lstsym(ntpsit,natms),nsym(ntpsit),unqsit(ntpsit))
      allocate(latm(ntpatm,natms),natm(ntpatm),unqatm(ntpatm))
      allocate(cpmdorder(natms),symorder(ntpsit),rootname(ind))
      allocate(numatoms(ntpmol))
c     real values.. not sure i need these any more
      read(nturd)alen,blen,clen,alpha,beta,gamma

c     read cell dimensions
      read(nturd)(cell(i),i=1,9)

c     read atomistic properties
      do i=1,natms
        read(nturd)xxx(i)
        read(nturd)yyy(i)
        read(nturd)zzz(i)
        read(nturd)atom(i)
        read(nturd)atmname(i)
        read(nturd)mass(i)
        read(nturd)molind(i)
        read(nturd)ltype(i)
      enddo

      read(nturd)(numatoms(i),i=1,ntpmol)
c     read symmetry related stuff
      read(nturd)lsym
 
      do i=1,ntpsit
        read(nturd)unqsit(i)
        read(nturd)nsym(i)
        read(nturd)(lstsym(i,j),j=1,nsym(i))
      enddo
      do i=1,ntpatm
        read(nturd)unqatm(i)
        read(nturd)natm(i)
        read(nturd)(latm(i,j),j=1,natm(i))
      enddo

      read(nturd)(rootname(i),i=1,ind)
      read(nturd)(title(i),i=1,80)
      read(nturd)(symorder(i),i=1,ntpsit) 

      read(nturd)(charges(i),i=1,ntpsit)
      read(nturd)lcpmdord
      if(lcpmdord)read(nturd)(cpmdorder(i),i=1,natms)
      if(lcpmdord)read(nturd)(printorder(i),i=1,natms)
   
      read(nturd)lvdw
      if(lvdw)then
        read(nturd)nvdw
        allocate(vdwatm(nvdw),eps(nvdw),sig(nvdw))
        read(nturd)(vdwatm(i),i=1,nvdw)
        read(nturd)(eps(i),i=1,nvdw)
        read(nturd)(sig(i),i=1,nvdw)
      endif 
      close(nturd)
      return
      end subroutine readturd
      subroutine images
     x  (natms,cell,xxx,yyy,zzz)
      
c***********************************************************************
c     
c     wl
c     2006/11/28 16:32:48
c     1.9
c     Exp
c     
c***********************************************************************
      
      implicit none

      integer natms,i
      real(8) cell,xxx,yyy,zzz,aaa,bbb,ccc,det,rt2,rt3,ssx
      real(8) ssy,ssz,ddd,xss,yss,zss,rcell

      dimension xxx(*),yyy(*),zzz(*)
      dimension cell(9),rcell(9)

      data rt2/1.41421356623d0/,rt3/1.7320508075d0/


c     parallelepiped boundary conditions
      call invert(cell,rcell,det)
         
      do i=1,natms
        ssx=(rcell(1)*xxx(i)+rcell(4)*yyy(i)+rcell(7)*zzz(i))
        ssy=(rcell(2)*xxx(i)+rcell(5)*yyy(i)+rcell(8)*zzz(i))
        ssz=(rcell(3)*xxx(i)+rcell(6)*yyy(i)+rcell(9)*zzz(i))
          
        xss=ssx-nint(ssx)
        yss=ssy-nint(ssy)
        zss=ssz-nint(ssz)
          
        xxx(i)=(cell(1)*xss+cell(4)*yss+cell(7)*zss)
        yyy(i)=(cell(2)*xss+cell(5)*yss+cell(8)*zss)
        zzz(i)=(cell(3)*xss+cell(6)*yss+cell(9)*zss)
          
      enddo
        
      return
      end subroutine images
      subroutine dcell(aaa,bbb)

c***********************************************************************
c     
c     dl_poly subroutine to calculate the dimensional properties of
c     a simulation cell specified by the input matrix aaa.
c     the results are returned in the array bbb, with :
c     
c     bbb(1 to 3) - lengths of cell vectors
c     bbb(4 to 6) - cosines of cell angles
c     bbb(7 to 9) - perpendicular cell widths
c     bbb(10)     - cell volume
c     
c     copyright daresbury laboratory 1992
c     author - w. smith         july 1992
c     
c     wl
c     2006/11/28 16:32:48
c     1.9
c     Exp
c     
c***********************************************************************

      implicit none

      real(8) aaa,bbb,axb1,axb2,axb3,bxc1,bxc2,bxc3,cxa1,cxa2,cxa3

      dimension aaa(9),bbb(10)

c     calculate lengths of cell vectors

      bbb(1)=sqrt(aaa(1)*aaa(1)+aaa(2)*aaa(2)+aaa(3)*aaa(3))
      bbb(2)=sqrt(aaa(4)*aaa(4)+aaa(5)*aaa(5)+aaa(6)*aaa(6))
      bbb(3)=sqrt(aaa(7)*aaa(7)+aaa(8)*aaa(8)+aaa(9)*aaa(9))

c     calculate cosines of cell angles

      bbb(4)=(aaa(1)*aaa(4)+aaa(2)*aaa(5)+aaa(3)*aaa(6))/(bbb(1)*bbb(2))
      bbb(5)=(aaa(1)*aaa(7)+aaa(2)*aaa(8)+aaa(3)*aaa(9))/(bbb(1)*bbb(3))
      bbb(6)=(aaa(4)*aaa(7)+aaa(5)*aaa(8)+aaa(6)*aaa(9))/(bbb(2)*bbb(3))

c     calculate vector products of cell vectors

      axb1=aaa(2)*aaa(6)-aaa(3)*aaa(5)
      axb2=aaa(3)*aaa(4)-aaa(1)*aaa(6)
      axb3=aaa(1)*aaa(5)-aaa(2)*aaa(4)
      bxc1=aaa(5)*aaa(9)-aaa(6)*aaa(8)
      bxc2=aaa(6)*aaa(7)-aaa(4)*aaa(9)
      bxc3=aaa(4)*aaa(8)-aaa(5)*aaa(7)
      cxa1=aaa(8)*aaa(3)-aaa(2)*aaa(9)
      cxa2=aaa(1)*aaa(9)-aaa(3)*aaa(7)
      cxa3=aaa(2)*aaa(7)-aaa(1)*aaa(8)

c     calculate volume of cell

      bbb(10)=abs(aaa(1)*bxc1+aaa(2)*bxc2+aaa(3)*bxc3)

c     calculate cell perpendicular widths

      bbb(7)=bbb(10)/sqrt(bxc1*bxc1+bxc2*bxc2+bxc3*bxc3)
      bbb(8)=bbb(10)/sqrt(cxa1*cxa1+cxa2*cxa2+cxa3*cxa3)
      bbb(9)=bbb(10)/sqrt(axb1*axb1+axb2*axb2+axb3*axb3)

      return
      end subroutine dcell
      subroutine invert(a,b,d)

c***********************************************************************
c     
c     dl_poly subroutine to invert a 3 * 3 matrix using cofactors
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith       april 1992
c     
c     wl
c     2006/11/28 16:32:48
c     1.9
c     Exp
c     
c***********************************************************************

      implicit none

      real(8) a,b,d,r

      dimension a(9),b(9)

c     calculate adjoint matrix
      b(1)=a(5)*a(9)-a(6)*a(8)
      b(2)=a(3)*a(8)-a(2)*a(9)
      b(3)=a(2)*a(6)-a(3)*a(5)
      b(4)=a(6)*a(7)-a(4)*a(9)
      b(5)=a(1)*a(9)-a(3)*a(7)
      b(6)=a(3)*a(4)-a(1)*a(6)
      b(7)=a(4)*a(8)-a(5)*a(7)
      b(8)=a(2)*a(7)-a(1)*a(8)
      b(9)=a(1)*a(5)-a(2)*a(4)

c     calculate determinant
      d=a(1)*b(1)+a(4)*b(2)+a(7)*b(3)
      r=0.d0
      if(abs(d).gt.0.d0)r=1.d0/d

c     complete inverse matrix
      b(1)=r*b(1)
      b(2)=r*b(2)
      b(3)=r*b(3)
      b(4)=r*b(4)
      b(5)=r*b(5)
      b(6)=r*b(6)
      b(7)=r*b(7)
      b(8)=r*b(8)
      b(9)=r*b(9)

      return
      end subroutine invert
      subroutine initscan(npdb,nmols,natms,ntpsit)
c*********************************************************************
c
c     subroutine reads pdb file for initial memory allocation
c
c*********************************************************************
      implicit none
      logical safe,done,molchk
      integer npdb,nmols,natms,ijunk,idum,imol,mol,ntpsit
      integer, dimension(5000) :: molcheck
      real(8) rjunk
      character*8 cjunk

      nmols=0
      data done/.false./
      natms=0
      do while(.not.done)
        call getrec(safe,npdb)
        call lowcase(record,lenrec)
        if(findstring('remark',record,idum))then
c          ignore
        elseif(findstring('atom',record,idum))then
          natms=natms+1
          ijunk=intstr(record,lenrec,idum)
          call getword(cjunk,record,3,idum)
          call getword(cjunk,record,3,idum)
          call getword(cjunk,record,3,idum)
          mol=dblstr(record,lenrec,idum)
c       check for unique molecules
          molchk=.true.
          do imol=1,nmols
           
            if(mol.eq.molcheck(imol))molchk=.false. 
          enddo
          if(molchk)then
            nmols=nmols+1
            molcheck(nmols)=mol
          endif

        elseif(findstring('end',record,idum))then
          done=.true.
        endif
        if(.not.safe)done=.true.
      enddo
      close(npdb)
      return
      end subroutine initscan
      subroutine readpdb(npdb,pdbfile,nmols,natms,lsym,ntpsit,
     & ntpatm,ntpmol,cell,alen,blen,clen,alpha,beta,gamma,title)
c*********************************************************************
c
c     subroutine reads pdb file and stores all the values in memory 
c
c*********************************************************************
      implicit none
      logical lsym,done,safe,sitchk,atchk,molchk
      integer npdb,nmols,natms,ntpsit,iatm,idum
      integer i,ijunk,mol,ntpatm,ntpmol,rmkcnt
      real(8) xcoord,ycoord,zcoord,rjunk
      real(8) alen,blen,clen,alpha,beta,gamma
      real(8), dimension(9) :: cell
      character*8 atname,sitnam,cjunk
      character*1 title(80)
      character*70 pdbfile

      done=.false.

      open(npdb,file=pdbfile,status='old')
      iatm=0
 
      rmkcnt=0
      ntpmol=0
      ntpsit=0
      ntpatm=0
      do while(.not.done)
        call getrec(safe,npdb)
        call lowcase(record,7)
        call strip(record,lenrec)
        if(findstring('remark',record,idum))then
          rmkcnt=rmkcnt+1
          if(rmkcnt.eq.1)then
            call getword(cjunk,record,6,lenrec)
            call copystring(record,title(1),80)
          endif
        elseif(findstring('cryst1',record,idum))then
          call getword(cjunk,record,6,lenrec)  
          alen=dblstr(record,lenrec,idum)
          blen=dblstr(record,lenrec,idum)
          clen=dblstr(record,lenrec,idum)
          alpha=dblstr(record,lenrec,idum)
          beta=dblstr(record,lenrec,idum)
          gamma=dblstr(record,lenrec,idum)
        elseif(findstring('atom',record,idum))then
          call getword(cjunk,record,6,lenrec) 
          iatm=iatm+1
          ijunk=intstr(record,lenrec,idum) 
          call getword(sitnam,record,3,lenrec) 
          call getword(cjunk,record,5,lenrec)

          mol=intstr(record,lenrec,idum)
c          numatoms(mol)=numatoms(mol)+1

          xcoord=dblstr(record,lenrec,idum)
          ycoord=dblstr(record,lenrec,idum)
          zcoord=dblstr(record,lenrec,idum)
 
          rjunk=dblstr(record,lenrec,idum)
          rjunk=dblstr(record,lenrec,idum)
          
          call getword(atname,record,3,lenrec)


          xxx(iatm)=xcoord
          yyy(iatm)=ycoord
          zzz(iatm)=zcoord

          atmname(iatm)=sitnam
          atom(iatm)=atname
          mass(iatm)=atmmass(atname)
c          molind(iatm)=mol

          
c        check for unique molecules
          molchk=.true.

          do i=1,ntpmol
            if(mol.eq.lstmol(i))then
              molchk=.false.
              numatoms(i)=numatoms(i)+1
              molind(iatm)=i
            endif
          enddo
          if(molchk)then
            ntpmol=ntpmol+1
            lstmol(ntpmol)=mol
            numatoms(ntpmol)=1
            molind(iatm)=ntpmol
          endif

      
c       check for unique atom types
      
          sitchk=.true.
          do i=1,ntpsit
    
            if(lsym.and.(sitnam(1:3).eq.unqsit(i)))then
              sitchk=.false.
              nsym(i)=nsym(i)+1
              lstsym(i,nsym(i))=iatm
              ltype(iatm)=i

            endif
          enddo

          if(sitchk)then
            ntpsit=ntpsit+1
            unqsit(ntpsit)=sitnam
            nsym(ntpsit)=1
            lstsym(ntpsit,1)=iatm
            ltype(iatm)=ntpsit
          endif 
       
c       check for unique atoms
          atchk=.true.
          do i=1,ntpatm
    
            if(atname(1:2).eq.unqatm(i))then
              atchk=.false.
              natm(i)=natm(i)+1
              latm(i,natm(i))=iatm

            endif
          enddo

          if(atchk)then
            ntpatm=ntpatm+1
            unqatm(ntpatm)=atname(1:2)
            natm(ntpatm)=1
            latm(ntpatm,1)=iatm
          endif 

        elseif(findstring('end',record,idum))then
          done=.true.
        endif
        if(.not.safe)done=.true.
      enddo

      if(lsym)call reorder(ntpmol,ntpsit,ntpatm,natms)

      close(npdb)
      return
      

      return
      end subroutine readpdb

      subroutine reorder(ntpmol,ntpsit,ntpatm,natms)
c*********************************************************************
c     re-order the atoms and their indices so that each molecule
c     is put together, and symmetry molecules are put one after the
c     other
c     * REALLY UGLY * - try to re-write this to make it simpler.
c*********************************************************************
      implicit none
      logical done,buffchk,ballschk,atchk,sitchk,molchk
      integer natms,i,j,k,m,n,tt,nbuff,iord,jatm,mol,ntpatm,catm,bufind
      integer atm,newatm,sityp,locmol,mm,mcnt,ntpmol,ntpsit,iatm,nrept
      
      integer, dimension(natms) :: buffer,neword,oldord,sitbuff
      integer, dimension(ntpmol) :: molly,atomcount
      integer, dimension(ntpmol,natms) :: symbuff

      real(8), dimension(natms) :: xbuf,ybuf,zbuf,masbuf
      integer, dimension(natms) :: molbuf
      character*3, dimension(natms) :: sitbuf
      character*2, dimension(natms) :: atombuf

      done=.false.

      do i=1,ntpmol
        molly(i)=0
        atomcount(i)=0
      enddo

      mcnt=0
c     iterate over symmetry types
      do i=1,ntpsit
c       iterate over all atoms
        do j=1,natms
c         check for like symmetry
          if(ltype(j).eq.i)then
            mol=molind(j)
            molchk=.true.
            do k=1,mcnt
              if(mol.eq.molly(k))then
                molchk=.false.
                bufind=k
              endif
            enddo
            if(molchk)then
              mcnt=mcnt+1
              molly(mcnt)=mol
              bufind=mcnt
            endif
            atomcount(bufind)=atomcount(bufind)+1
            symbuff(bufind,atomcount(bufind))=j
          endif
        enddo
      enddo

c     now we should have mcnt equal to ntpmol, atomcount should have
c     the same number of atoms as the equvalent numatoms (different 
c     order), and symbuff should contain the re-ordered atoms

      if(mcnt.ne.ntpmol)then
        write(*,"('ERROR - when re-ordering the atoms the number of 
     &molecules was miscounted.  Counted ',i3,' should be ',i3)")
     &mcnt,ntpmol
      endif
c     convert to 1d
      iatm=0
      do i=1,mcnt
        do j=1,atomcount(i)
          iatm=iatm+1
          jatm=symbuff(i,j)
          oldord(iatm)=jatm
          neword(jatm)=iatm
        enddo  

      enddo

c     re-order the arrays....
      do i=1,natms
        newatm=oldord(i)
  
        xbuf(i)=xxx(newatm)
        ybuf(i)=yyy(newatm)
        zbuf(i)=zzz(newatm)
        masbuf(i)=mass(newatm)
        sitbuf(i)=atmname(newatm)
        atombuf(i)=atom(newatm)
        molbuf(i)=molind(newatm)

      enddo
c     re-number the molecular indices

      xxx=xbuf
      yyy=ybuf
      zzz=zbuf
      mass=masbuf
      atmname=sitbuf
      atom=atombuf
      molind=molbuf
c     check for unique molecules

      do i=1,ntpsit
        nsym(i)=0
      enddo
      do i=1,ntpmol
        numatoms(i)=0
      enddo
      do i=1,ntpatm
        natm(i)=0
      enddo
      ntpmol=0
      ntpatm=0
      ntpsit=0
      do j=1,natms
        molchk=.true.
        do i=1,ntpmol
          if(molind(j).eq.lstmol(i))then
            molchk=.false.
            numatoms(i)=numatoms(i)+1
            molind(j)=i
          endif
        enddo
        if(molchk)then
          ntpmol=ntpmol+1
          lstmol(ntpmol)=molind(j)
          numatoms(ntpmol)=1
          molind(j)=ntpmol
        endif

     
c     check for unique atom types 
        sitchk=.true.
        do i=1,ntpsit
     
          if(atmname(j).eq.unqsit(i))then
            sitchk=.false.
            nsym(i)=nsym(i)+1
            lstsym(i,nsym(i))=j
            ltype(j)=i

          endif
        enddo

        if(sitchk)then
          ntpsit=ntpsit+1
          unqsit(ntpsit)=atmname(j)
          nsym(ntpsit)=1
          lstsym(ntpsit,1)=j
          ltype(j)=ntpsit
        endif 
c     check for unique atoms
        atchk=.true.
        do i=1,ntpatm
     
          if(atom(j).eq.unqatm(i))then
            atchk=.false.
            natm(i)=natm(i)+1
            latm(i,natm(i))=j

          endif
        enddo

        if(atchk)then
          ntpatm=ntpatm+1
          unqatm(ntpatm)=atom(j)
          natm(ntpatm)=1
          latm(ntpatm,1)=j
        endif 

      enddo 
c      print *, atmname(2)
c      do i=1,iord
c        print *, oldord(i),neword(i),oldord(neword(i)) 
c      enddo
      return
      end subroutine reorder

      subroutine argparse(pdbfile,init,lcpmd,cpmdread,reptwrite,
     &reptread,geom,elpot,lvdw,fldwrite,na,nb,nc,cnfgwrite,lpdb
     &,ljfile,lsym)
c*********************************************************************
c
c     subroutine reads in all command line arguments and returns
c     logical arguments necessary to proceed with stuff
c
c*********************************************************************
      implicit none
      logical init,lcpmd,geom,elpot,cpmdread,reptwrite,reptread
      logical lvdw,fldwrite,cnfgwrite,lpdb,lsym
      integer nargs,i,n,j,idum,na,nb,nc
      character*70 arg_str,arg_two,pdbfile,ljfile
      
      nargs=iargc()

      if(nargs.eq.0)then
        write(*,"('ERROR - no arguments defined on command line')")
        stop
      endif
      j=0
      do n=1,nargs
        j=j+1
        call getarg(j,arg_str)
c       convert argument to array
        do i=1,70
          record(i)=arg_str(i:i)
        enddo
        call lowcase(record,lenrec)


        if(findstring('-f',record,idum))then
          j=j+1
          nargs=nargs-1
          call getarg(j,pdbfile)
          init=.true.
        elseif(findstring('-vdw',record,idum))then
          j=j+1
          lvdw=.true.
          call getarg(j,ljfile)
        elseif(findstring('-pdb',record,idum))then
          lpdb=.true.
        elseif(findstring('-nosym',record,idum))then
          lsym=.false.
        elseif(findstring('-dlp',record,idum))then
          cnfgwrite=.true.
          fldwrite=.true.
          j=j+1
          call getarg(j,arg_str)
          do i=1,70
            record(i)=arg_str(i:i)
          enddo
          na=intstr(record,70,idum)
          j=j+1
          call getarg(j,arg_str)
          do i=1,70
            record(i)=arg_str(i:i)
          enddo
          nb=intstr(record,70,idum)
          j=j+1
          call getarg(j,arg_str)
          do i=1,70
            record(i)=arg_str(i:i)
          enddo
          nc=intstr(record,70,idum)
          if(na.eq.0)na=1
          if(nb.eq.0)nb=1
          if(nc.eq.0)nc=1

        elseif(findstring('-debug',record,idum))then
        elseif(findstring('-cpmd',record,idum))then
          j=j+1
          
          call getarg(j,arg_two)
          if(arg_two.eq.'geom')then
            lcpmd=.true.
            geom=.true.
          elseif(arg_two.eq.'elpot')then
            lcpmd=.true.
            elpot=.true.
          elseif(arg_two.eq.'read')then
            cpmdread=.true.
          else
            write(*,"('ERROR - unrecognized cpmd option :',a70)")
     &arg_two
          endif
        elseif(findstring('-rept',record,idum))then
          j=j+1
         
          call getarg(j,arg_two)
          if(arg_two.eq.'read')then
            reptread=.true.
          elseif(arg_two.eq.'write')then
            reptwrite=.true.
          endif 
        elseif(arg_str.eq.' ')then
        else
          write(*,"('ERROR - cannot recognize argument :',a70)")arg_str
          stop
        endif
      enddo

      return
      end subroutine argparse

      subroutine getcell(cell,alen,blen,clen,alpha,beta,gam) 
c*********************************************************************
c     
c     subroutine to convert to cell vectors 
c     from alen,blen,clen
c     alpha = angle between b and c
c     beta  = angle between a and c
c     gamma = angle between a and b
c*********************************************************************
      implicit none
      integer i
      real(8), parameter :: deg2rad=0.0174532925d0
      real(8), parameter :: rad2deg=57.2957795d0
      real(8), parameter :: pi=3.141592653589793d0
      real(8) alen,blen,clen,alpha,beta,gam
      real(8) pr
      real(8), dimension(9) :: cell
      real(8), dimension(3) :: avec,bvec,cvec,v1,v2,v3

      alpha=alpha*deg2rad
      beta=beta*deg2rad
      gam=gam*deg2rad
c     the a axis is oriented to the xaxis
c     start with orthogonal basis
      v1=(/1.d0,0.d0,0.d0/)

      v2=(/0.d0,1.d0,0.d0/)

      v3=(/0.d0,0.d0,1.d0/)

      avec=v1

      bvec=dsin(gam)*v2+dcos(gam)*v1

      pr=(dcos(alpha)-dcos(gam)*dcos(beta))/dsin(gam)

      cvec=v3*sqrt(1-dcos(beta)**2-pr**2)
     &+dcos(beta)*v1+pr*v2
c     write to cell
      cell(1:3)=avec*alen
      cell(4:6)=bvec*blen
      cell(7:9)=cvec*clen
c     fix for rounding
      do i=1,9
        cell(i)=float(nint(cell(i)*1.d6))/1.d6
      enddo

      return
      end subroutine getcell
      subroutine getrec(safe,ifile)

c*********************************************************************
c     
c     dl_poly subroutine to read a character string on one node
c     
c*********************************************************************

      implicit none
      
      logical safe

      character*150 line
      integer ifile,i
      
      safe=.true.
      
      read(ifile,'(a150)',end=100)line
 
      do i=1,lenrec

        record(i)=line(i:i)
         
      enddo
             
      return
        
  100 safe=.false.
                  
      end subroutine getrec

      real(8) function atmmass(i) 
c******************************************************
c     generate an atom label from an atomic number 
c
c******************************************************
      implicit none
      character*2 i
      if (i.eq.'H ')then
        atmmass=1.00794
      elseif(i.eq.'He')then
        atmmass=4.002602
      elseif(i.eq.'Li')then
        atmmass=6.941
      elseif(i.eq.'Be')then
        atmmass=9.012182
      elseif(i.eq.'B ')then
        atmmass=10.811
      elseif(i.eq.'C ')then
        atmmass=12.0107
      elseif(i.eq.'N ')then
        atmmass=14.00674
      elseif(i.eq.'O ')then
        atmmass=15.9994
      elseif(i.eq.'F ')then
        atmmass=18.9984032
      elseif(i.eq.'Ne')then
        atmmass=20.1797
      elseif(i.eq.'Na')then
        atmmass=22.989770
      elseif(i.eq.'Mg')then
        atmmass=24.3050
      elseif(i.eq.'Al')then
        atmmass=26.981538
      elseif(i.eq.'Si')then
        atmmass=28.0855
      elseif(i.eq.'P ')then
        atmmass=30.973761
      elseif(i.eq.'S ')then
        atmmass=32.066
      elseif(i.eq.'Cl')then
        atmmass=35.4527
      elseif(i.eq.'Ar')then
        atmmass=39.948
      elseif(i.eq.'K ')then
        atmmass=39.0983
      elseif(i.eq.'Ca')then
        atmmass=40.078
      elseif(i.eq.'Sc')then
        atmmass=44.955910
      elseif(i.eq.'Ti')then
        atmmass=47.867
      elseif(i.eq.'V ')then
        atmmass=50.9415
      elseif(i.eq.'Cr')then
        atmmass=51.9961
      elseif(i.eq.'Mn')then
        atmmass=54.938049
      elseif(i.eq.'Fe')then
        atmmass=55.845
      elseif(i.eq.'Co')then
        atmmass=58.9332
      elseif(i.eq.'Ni')then
        atmmass=58.6934
      elseif(i.eq.'Cu')then
        atmmass=63.546
      elseif(i.eq.'Zn')then
        atmmass=65.39
      elseif(i.eq.'Ga')then
        atmmass=69.723
      elseif(i.eq.'Ge')then
        atmmass=72.61
      elseif(i.eq.'As')then
        atmmass=74.9216
      elseif(i.eq.'Se')then
        atmmass=78.96
      elseif(i.eq.'Br')then
        atmmass=79.904
      elseif(i.eq.'Kr')then
        atmmass=83.80
      elseif(i.eq.'Y ')then
        atmmass=88.90585
      elseif(i.eq.'Zr')then
        atmmass=91.224
      elseif(i.eq.'Nb')then
        atmmass=92.90638
      elseif(i.eq.'Mo')then
        atmmass=95.94
      elseif(i.eq.'Ru')then
        atmmass=101.07
      elseif(i.eq.'Rh')then
        atmmass=102.90550
      elseif(i.eq.'Pd')then
        atmmass=106.42
      elseif(i.eq.'Ag')then
        atmmass=107.8682
      elseif(i.eq.'Cd')then
        atmmass=112.411
      elseif(i.eq.'In')then
        atmmass=114.818
      elseif(i.eq.'Sn')then
        atmmass=118.710
      elseif(i.eq.'Sb')then
        atmmass=121.760
      elseif(i.eq.'Te')then
        atmmass=127.760
      elseif(i.eq.'I ')then
        atmmass=126.90447
      elseif(i.eq.'Xe')then
        atmmass=131.29
      elseif(i.eq.'Hf')then
        atmmass=178.49
      elseif(i.eq.'Ta')then
        atmmass=180.9479
      elseif(i.eq.'W ')then
        atmmass=183.84
      elseif(i.eq.'Re')then
        atmmass=186.207
      elseif(i.eq.'Os')then
        atmmass=190.23
      elseif(i.eq.'Ir')then
        atmmass=192.217
      elseif(i.eq.'Pt')then
        atmmass=195.078
      elseif(i.eq.'Au')then
        atmmass=196.96655
      elseif(i.eq.'Hg')then
        atmmass=200.59
      elseif(i.eq.'Tl')then
        atmmass=204.3833
      elseif(i.eq.'Pb')then
        atmmass=207.2
      endif
      return
      end function atmmass
      integer function intstr(word,len,lst)

c***********************************************************************
c     
c     dl_poly function for extracting integers from a 
c     character string
c     
c     copyright - daresbury laboratory 1994
c     author    - w. smith may 1994.
c     
c     parameters:
c     word   - input character string
c     len    - working length of character string
c     lst    - location of space character at end of
c     integer string
c     
c     wl
c     2007/07/31 10:02:58
c     1.3
c     Exp
c     
c***********************************************************************
      
      implicit none

      logical flag,count,final
      character*1 n,word,ksn
      integer lst,len,j,isn

      dimension n(0:9),word(len)
      data n/'0','1','2','3','4','5','6','7','8','9'/

      isn=1
      lst=0
      ksn='+'
      intstr=0
      flag=.false.
      final=.false.
      count=.false.
      
      do while(lst.lt.len.and.(.not.final))

        lst=lst+1
        flag=.false.

        do j=0,9
          
          if(n(j).eq.word(lst))then
            
            intstr=10*intstr+j
            count=.true.
            flag=.true.
            
          endif
          
        enddo

        if(count.and.(.not.flag))final=.true.
        if(flag.and.ksn.eq.'-')isn=-1
        ksn=word(lst)

      enddo

      intstr=isn*intstr

      do j=lst,len
        word(j-lst+1)=word(j)
      enddo
      do j=len-lst+2,len
        word(j)=' '
      enddo

      return
      end function intstr

      real(8) function dblstr(word,len,lst)

c***********************************************************************
c     
c     dl_poly function for extracting double precisions from a 
c     character string. 
c     modified from dl_poly function intstr
c     
c     copyright - daresbury laboratory 1994
c     author    - w. smith may 1994.
c     modified  - t. forester april 1994
c     
c     parameters:
c     word   - input character string
c     len    - working length of character string
c     lst    - location of space character at end of
c     double precision string
c     
c     wl
c     2007/07/31 10:02:58
c     1.3
c     Exp
c     
c***********************************************************************
      
      implicit none
      
      character*1 n,word,ksn,dot,d,e
      logical flag,ldot,start,final
      integer len,lst,iexp,idum,i,j,fail
      real(8) sn,ten,one
      dimension n(0:9),word(len)
      character*1, allocatable :: work(:)

      data n/'0','1','2','3','4','5','6','7','8','9'/
      data dot/'.'/
      data d/'d'/
      data e/'e'/
      
      allocate(work(len),stat=fail)

      lst=0
      sn=1.d0
      ksn='+'
      ten=10.d0
      one=1.d0
      
      dblstr=0.d0
      iexp=0
      idum=0
      start=.false.
      ldot=.false.
      final=.false.

      do while(lst.lt.len.and.(.not.final))
        
        lst=lst+1
        flag=.false.
        
        do j=0,9
          
          if(n(j).eq.word(lst))then
            
            dblstr=ten*dblstr+one*dble(j)
            flag=.true.
            start=.true.
            
          endif
          
        enddo
        
        if(dot.eq.word(lst))then
          
          flag=.true.
          ten=1.d0
          ldot=.true.
          start=.true.
          
        endif

        if(flag.and.ksn.eq.'-') sn=-1.d0
        if(ldot) one=one/10.d0
        ksn=word(lst)
        if(ksn.eq."D")ksn="d"
        if(ksn.eq."E")ksn="e"
        
        if(start)then
          
          if(d.eq.ksn.or.e.eq.ksn)then
            
            do i=1,len-lst
              work(i)=word(i+lst)
            enddo
            iexp=intstr(work,len-lst,idum)
            final=.true.

          endif

          if(.not.flag)final=.true.
          
        endif
        
      enddo
      
      dblstr=sn*dblstr*(10.d0**iexp)
      lst=lst+idum
      
      do j=lst,len
        word(j-lst+1)=word(j)
      enddo
      do j=len-lst+2,len
        word(j)=' '
      enddo

      deallocate(work,stat=idum)

      return
      end function dblstr

      subroutine strip(string,imax)

c***********************************************************************
c     
c     DL_POLY routine to strip blanks from start of a string
c     maximum length is 255 characters
c     
c     copyright daresbury laboratory 1993
c     author   t.forester       july 1993
c     
c     wl
c     2007/07/31 10:02:58
c     1.3
c     Exp
c     
c***********************************************************************

      implicit none

      integer i,imax,j
      character*1 string(imax)
      
      do i=1,imax

        if(string(1).eq.' ')then

          do j=1,imax-1

            string(j)=string(j+1)

          enddo

          string(imax)=' '

        endif

      enddo

      return
      end subroutine strip

      subroutine lowcase(string,length)

c***********************************************************************
c     
c     DL_POLY routine to lowercase a string of up to 255 characters.
c     Transportable to non-ASCII machines
c     
c     copyright daresbury laboratory 1993
c     author    t. forester     july 1993
c     
c     wl
c     2007/07/31 10:02:58
c     1.3
c     Exp
c     
c***********************************************************************

      implicit none

      character*1 string(*)
      character*1 letter
      integer i,length

      do i=1,min(255,length)

        letter=string(i)

        if(letter.eq.'A')then
          letter='a'
        else if(letter.eq.'B')then
          letter='b'
        else if(letter.eq.'C')then
          letter='c'
        else if(letter.eq.'D')then
          letter='d'
        else if(letter.eq.'E')then
          letter='e'
        else if(letter.eq.'F')then
          letter='f'
        else if(letter.eq.'G')then
          letter='g'
        else if(letter.eq.'H')then
          letter='h'
        else if(letter.eq.'I')then
          letter='i'
        else if(letter.eq.'J')then
          letter='j'
        else if(letter.eq.'K')then
          letter='k'
        else if(letter.eq.'L')then
          letter='l'
        else if(letter.eq.'M')then
          letter='m'
        else if(letter.eq.'N')then
          letter='n'
        else if(letter.eq.'O')then
          letter='o'
        else if(letter.eq.'P')then
          letter='p'
        else if(letter.eq.'Q')then
          letter='q'
        else if(letter.eq.'R')then
          letter='r'
        else if(letter.eq.'S')then
          letter='s'
        else if(letter.eq.'T')then
          letter='t'
        else if(letter.eq.'U')then
          letter='u'
        else if(letter.eq.'V')then
          letter='v'
        else if(letter.eq.'W')then
          letter='w'
        else if(letter.eq.'X')then
          letter='x'
        else if(letter.eq.'Y')then
          letter='y'
        else if(letter.eq.'Z')then
          letter='z'
        endif

        string(i)=letter

      enddo

      return
      end subroutine lowcase

      subroutine copystring(oldstr,newstr,length)

c***********************************************************************
c     
c     DL_POLY routine to copy one string into another
c     
c     copyright daresbury laboratory
c     author    w. smith    jan 2004
c     
c     wl
c     2007/07/31 10:02:58
c     1.3
c     Exp
c     
c***********************************************************************

      implicit none

      character*1 newstr(*),oldstr(*)
      integer i,length

      do i=1,length

        newstr(i)=oldstr(i)

      enddo

      return
      end subroutine copystring

      logical function findstring(seek,string,here)

c***********************************************************************
c     
c     DL_POLY routine to find an explicit string in an input record
c     note: variable `seek' is a character string while variable
c    `string' is a character*1 array i.e. code is application specific
c
c     copyright daresbury laboratory
c     author    w.smith   jan   2004
c     
c     wl
c     2007/07/31 10:02:58
c     1.3
c     Exp
c     
c***********************************************************************

      implicit none

      integer i,n,m,here
      character*(*) seek
      character*1 string(lenrec)

      m=lenrec
      n=len(seek)
      findstring=.false.

      here=0
      do while(here.le.m-n.and.(.not.findstring))

        findstring=.true.

        do i=1,n
          if(seek(i:i).ne.string(here+i))findstring=.false.
        enddo

        here=here+1

      enddo

      return
      end function findstring

      subroutine getword(word,string,len1,len2)

c***********************************************************************
c     
c     DL_POLY routine to fetch an 8 character word from a string
c     while ignoring leading blanks
c
c     copyright daresbury laboratory
c     author   w.smith jan 2004
c     
c     wl
c     2007/07/31 10:02:58
c     1.3
c     Exp
c     
c***********************************************************************

      implicit none

      logical final
      character*8 word
      integer len1,len2,i,j,k
      character*1 wrdseq(len1),string(len2)
c      character*1 word(len1),string(len2)
      do i=1,len1
        wrdseq(i)=' '
c         word(i)=' '
      enddo

      i=0
      k=0
      final=.false.
      
      do while(.not.final.and.i.lt.len2)
        
        i=i+1
        
        if(string(1).eq.' ')then
          
          if(k.gt.0)final=.true.
          
        else
          
          k=k+1
          wrdseq(k)=string(1)
c          word(k)=string(1)
          if(k.eq.len1)final=.true.

        endif
        
        do j=1,len2-1
          
          string(j)=string(j+1)
          
        enddo
        
        string(len2)=' '
          
      enddo
      
      word=mkwd8(wrdseq)

      return
      end subroutine getword

      character*8 function mkwd8(string)

c***********************************************************************
c     
c     DL_POLY routine to make an 8 character word from a string
c
c     copyright daresbury laboratory
c     author   w.smith nov 2006
c     
c     wl
c     2007/07/31 10:02:58
c     1.3
c     Exp
c     
c***********************************************************************

      implicit none

      integer i
      character*1 string(*)
      
      do i=1,8
         mkwd8(i:i)=string(i)
      enddo
      
      return
      end function mkwd8
      
      end program 
