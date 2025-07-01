        program Traff_all  

c     Model purpose
c     Called FLEE (Forecasting Laboratory for Exploring the Evacuation-system), this agent-based model simulates the evacuation response to hurricane forecasts. 
c     In particular, it serves as a virtual laboratory to explore the dynamics of the hurricane forecast-warning-evacuation system 

c     Model inputs
c     Forecasts- Every 6 hours in the model, this program will ingest forecast files (excel inputs e.g., ADV26), located in the same folder as this .f file.
c     The files contain the wind, surge, and rain risk across the model, presented on a 0-10 scale. Collectively, they simulate an evolving hurricane forecast scenario. 
c     Risk files are seperately generated from archived forecasts from the National Hurricane Center including the hurricane's track, intensity, and the cone of uncertainty. 

c     Built environment - A second input is the npop.txt file. This has information on the "people" side of the model, such as how many people reside in the model and where. Other info includes 
c     the amount of time it takes to evacuate each grid cell (clerance times), inundation potential (based on NHC storm surge data), and the average characteristics of people in a grid cell
c     (based on social vulnerability data form the CDC). These variables are static, meaning they are only defined once.  

c     The virtual world
c     A NXW x NYH grid of tile (8x20) representing the N-S axis of the Florida peninsula. The population, based on numbers in the
c     npop file, are geographically distributed in such a way to simulate the actual population of Florida. For example, the population centers of Miami and Tampa are
c     distinguishable in the model. To create this virtual world, geospatial population and social datasets were projected onto the agent-based model grid using GIS. 
c     While idealized, it makes the FLEE's virtual world "Florida-like" 

c     Agents
c     There are three primary types of agents in FLEE. 

c     (1) Emergency manager agents. For each grid cell, there is one emergency manager agent who make a  
c     decision about whether or not to issue evacuation order for their grid cell (and when). That decision is revisited every hour, based on the forecast risk at its grid cell, 
c     the estimated time of arrival of the hurricane, and their clearance times (estimated amount of time it takes to evacuate that grid cell). 

c     (2) Household agents. Each household (units of 4 people) has a unique set of characteristics known to influence evacuation decision-making. These include (a) forecast risk, 
c     (b) evacuation orders, (c) mobile home ownership, (d) car ownership, and (e) the age of the household. For c-e, household characteristics are stochastically assigned (only once, in the subroutine)
c     in this code based on the average grid cell-level characteristcs in the npop file. 
c     
c     Household agents make decisions about if they should evacuate, a decision they revisit every 30 minutes. 
c     That decision process is located in the subroutine at the bottom of the code. 
c     If a household decides to evacuate, they also decide where to go (a decision made once and immediately afterward). They also decide when to leave, a decision that's influenced by the storms arrival time. 
c     These decision making processes are rooted in empirical social science studies. 

c     (3) Vehicle agents. If a household decides to evacuate, they become a vehicle agent. These agents are initially placed onto an open spot on a nearby road. 
c     Vehicles accelerate when they can, maintain speed at the speed limit or when behind another car, and slow down if they must.  

c     Traffic model
c     The evacuation traffic model dumps out the state at every minute of
c     integration (50 time steps using an increment of 1.2 sec). This high-res aspect of the model is needed to simulate vehicle-vehicle interactions
c     which can generate traffic jams. 

c     Built environment 
c     As part of the Florida-like virtual world, idealized road networks are placed in between grid cells, some moving E-W, others moving S-N. 

c ... N/S Interstates:
c     At the junction of X=8 (east Florida) is a N/S -running Interstate (5 lanes) extending from Y=1-NY where 1 is south. This is like I-95 running along FL's east coast.
c     At the junction of X=1 (west Florida) is an identential N/S running, 5 lane interstate extending to NY=20. This is like I-75 on Florida's west coast. 

c ... E/W Interstates
c     At the base of the grid Y=1 is an E-W running Interstate (3 lanes each way). These connect the two outer interstates/ 
c.    At the middle of the grid Y=10 is another E-W running Interstate (3 lanes each way). This connects Tampa with Orlando. 

c ... E/W Highways 
c.    In between each grid cell are E-W running higways (2-lanes) that  joins the traffic from the W-E routes
c     to the N route.

c     Each highway consists of NCELL=3500 road-lengths (road length = 7.5m so NCELL = 26.25 km)
c     This means the grid is about 110 km wide and 262 km tall -- roughly FL panhandle-like.  

c     Setting initial parameters and the sizes of arrays
c     We set the interstate speed limit (nVMX) at 5 road-lengths per time step (1.2 s) which is about 70 mph. 
c     We set the highway speed limit (nVMX2) to 3 road-lengths per time step (1.2s) which is about 45 mph
c     We place random accidents (RA) on the major interstates every 2.7 hours or so. 
c     About 35% of cars will exhibit some "chaotic" or random movements (RP) which can trigger abrubt slowdowns/jams

       PARAMETER (NCELLX=1750,NCELLY=35000,NX=8,NY=20,nVMX=5,nVMX2=3,
     *   RLC=7.5,DT=1.2,RP=0.05,RA=0.00005)

c ... for the right (east) interstate....5 indicates the number of possible lanes     
       integer*4 TTime,iseed
       integer car(NCELLY,5),vel(NCELLY,5),lag(NCELLY,5)
       integer car2(NCELLY,5),vel2(NCELLY,5),ident2(NCELLY,5)
       integer carp2(NCELLY,5),carp(NCELLY,5),rt2(NCELLY,5)
       integer twoe(NCELLY,5),exit2(NCELLY,5),exite(NCELLY,5)
       integer two2e(NCELLY,5),ident(NCELLY,5)
c ... for the left (west) interstate
       integer carX(NCELLY,5),velX(NCELLY,5),lagX(NCELLY,5)
       integer car2X(NCELLY,5),vel2X(NCELLY,5),ident2X(NCELLY,5)
       integer carp2X(NCELLY,5),identX(NCELLY,5),carpX(NCELLY,5)
       integer one2w(NCELLY,5),two2w(NCELLY,5),exit2w(NCELLY,5)
       integer onew(NCELLY,5),twow(NCELLY,5),exitw(NCELLY,5)
c ... for the bottom interstate (e running)... 3 indicates the number of possible lanes
       integer carse(NX,NY,NCELLX,3),velse(NX,NY,NCELLX,3)
       integer lagse(NX,NY,NCELLX,3),identse(NX,NY,NCELLX,3)
       integer car2se(NX,NY,NCELLX,3),vel2se(NX,NY,NCELLX,3)
       integer ident2se(NX,NY,NCELLX,3)
       integer carp2se(NX,NY,NCELLX,3),carpse(NX,NY,NCELLX,3)
       integer onese(NX,NY,NCELLX,3),one2se(NX,NY,NCELLX,3)
       integer twose(NX,NY,NCELLX,3),two2se(NX,NY,NCELLX,3)
       integer exitse(NX,NY,NCELLX,3),exit2se(NX,NY,NCELLX,3)
c ... for the bottom interstate (w running)
       integer carsw(NX,NY,NCELLX,3),velsw(NX,NY,NCELLX,3)
       integer lagsw(NX,NY,NCELLX,3),identsw(NX,NY,NCELLX,3)
       integer car2sw(NX,NY,NCELLX,3),vel2sw(NX,NY,NCELLX,3)
       integer ident2sw(NX,NY,NCELLX,3)
       integer carp2sw(NX,NY,NCELLX,3),carpsw(NX,NY,NCELLX,3)      
       integer onesw(NX,NY,NCELLX,3),one2sw(NX,NY,NCELLX,3)  
       integer twosw(NX,NY,NCELLX,3),two2sw(NX,NY,NCELLX,3)       
       integer exitsw(NX,NY,NCELLX,3),exit2sw(NX,NY,NCELLX,3)
c ... for the middle interstate (w running)
       integer carmw(NX,NY,NCELLX,3),velmw(NX,NY,NCELLX,3)
       integer lagmw(NX,NY,NCELLX,3),identmw(NX,NY,NCELLX,3)
       integer car2mw(NX,NY,NCELLX,3),vel2mw(NX,NY,NCELLX,3)
       integer ident2mw(NX,NY,NCELLX,3)
       integer carp2mw(NX,NY,NCELLX,3),carpmw(NX,NY,NCELLX,3)      
       integer onemw(NX,NY,NCELLX,3),one2mw(NX,NY,NCELLX,3)  
       integer twomw(NX,NY,NCELLX,3),two2mw(NX,NY,NCELLX,3)       
       integer exitmw(NX,NY,NCELLX,3),exit2mw(NX,NY,NCELLX,3) 
c ... for the middle interstate (e running)
       integer carme(NX,NY,NCELLX,3),velme(NX,NY,NCELLX,3)
       integer lagme(NX,NY,NCELLX,3),identme(NX,NY,NCELLX,3)
       integer car2me(NX,NY,NCELLX,3),vel2me(NX,NY,NCELLX,3)
       integer ident2me(NX,NY,NCELLX,3)
       integer carp2me(NX,NY,NCELLX,3),carpme(NX,NY,NCELLX,3)      
       integer oneme(NX,NY,NCELLX,3),one2me(NX,NY,NCELLX,3)  
       integer twome(NX,NY,NCELLX,3),two2me(NX,NY,NCELLX,3)       
       integer exitme(NX,NY,NCELLX,3),exit2me(NX,NY,NCELLX,3)        
c ... for the highways ... 2 indicates the number of possible lanes
       integer care(NX,NY,NCELLX,2),vele(NX,NY,NCELLX,2)
       integer lage(NX,NY,NCELLX,2),idente(NX,NY,NCELLX,2)
       integer car2e(NX,NY,NCELLX,2),vel2e(NX,NY,NCELLX,2)
       integer ident2e(NX,NY,NCELLX,2),npop(NX,NY),npop2(NX,NY)
       integer carp2e(NX,NY,NCELLX,2),carpe(NX,NY,NCELLX,2)
       integer oneh(NX,NY,NCELLX,2),one2h(NX,NY,NCELLX,2)
       integer twoh(NX,NY,NCELLX,2),two2h(NX,NY,NCELLX,2)
       integer exith(NX,NY,NCELLX,2),exit2h(NX,NY,NCELLX,2)

c... for the agent characteristics 
c    We have N agents corresponding to total population in each grid cell. Largest tile pop is 
c    1,181,420 (TB) and overall population is 19,755,788. 
c    We can combine populations into households of 4 that behaves collectively.
c    Then the array size needed for each grid cell would be about 295,355.
c    Across all grid cells, that's 4,938,947. 
       integer orisk(NX,NY),wrisk(NX,NY),srisk(NX,NY),rrisk(NX,NY)
       integer storm(NX,NY,295355)
       integer evacthres(NX,NY),evac(NX,NY)
       integer socio(NX,NY),age(NX,NY),nocar(NX,NY),mobl(NX,NY)
       integer socios(NX,NY,295355),ages(NX,NY,295355)
       integer mobls(NX,NY,295355),agewt(NX,NY,295355)       
       integer ewt(NX,NY,295355)
       integer soph(4938947),pevac(NX,NY,295355)
       integer threatC(4938947),gender(4938947),priorEv(4938947)
       integer priorEx(4938947),priorCM(4938947),mevac(NX,NY,295355)
       INTEGER failatt(NX,NY,295355),dthresh(NX,NY,295355)
       INTEGER decisiont(NX,NY,295355),dest(NX,21)
       INTEGER destside(NX,NY,295355),onerte(NX,NY,295355)
       INTEGER tworte(NX,NY,295355),exitrow(NX,NY,295355)
       integer wwt(NX,NY,295355),swt(NX,NY,295355),rwt(NX,NY,295355)          
       integer barr(NX,NY,295355),awt(NX,NY,295355),mobwt(NX,NY,295355) 
       integer strmwt(NX,NY,295355),ct12(NX,NY),dect(NX,NY,295355)  
c    Defining arrays to count/track evacuation behaviors in each cell 	   
       INTEGER ct0(NX,NY),ct1(NX,NY),ct2(NX,NY),ct3(NX,NY),ct4(NX,NY)
       INTEGER ct5(NX,NY),ct6(NX,NY),maxspots(NX,21),spt(NX,21)
       integer ct7(NX,NY),ct8(NX,NY),ct9(NX,NY),ct10(NX,NY),ct11(NX,NY)
       integer ct99(NX,NY)
       integer ctimey(NX,NY),ctimeo(NX,NY),ctimer(NX,NY),tstoa(NX,NY)
       integer eotime(NX,NY),parteos(NX,NY),eoticks(NX,NY)
       
       iseed=Time()
       call srand(iseed)
       ivac=1

c ... open text file with grid population and household characterstics averaged across the grid. This is done once.    
c.     npop=number of people in each grid cell. 
c.     evacthres notes which grid cells are "eligeble" for evacuation orders (0=inland)
c.     socio, age, nocar, mobl are placed on a 1-5 scale based on social vulnerability datasets across the state
c.     ctimey, ctimeo, ctimer denote three levels of clearance times (low, middle, high) which are adjusted based on surge risk 
c.     parteo notes the percentage of a grid cell that is in an inundation or evacuation zone 

       open(unit=1,file='npop.txt',status='old')
       do while(.true.)
         read(1,*,end=50) mx,my,npop(mx,my),evacthres(mx,my),
     *       socio(mx,my),age(mx,my),nocar(mx,my),mobl(mx,my),
     *       ctimey(mx,my),ctimeo(mx,my),ctimer(mx,my),
     *       parteos(mx,my)
         npop2(mx,my)=npop(mx,my)
       end do      
 50    close(1)
 
 c.     Restructuring the age distribution of grid cells such as middle ages are more likely to evacuate. Could easily be done in the input files but 
 c.     quickly did them here. 
        do kx=1,8
          do ky=1,20       
             evac(kx,ky)=0
			 eoticks(kx,ky)=0
             if (age(kx,ky).eq.1) then
                age(kx,ky)=1
             else if (age(kx,ky).eq.2) then
                age(kx,ky)=3
             else if (age(kx,ky).eq.3) then
                age(kx,ky)=5
             else if (age(kx,ky).eq.4) then
                age(kx,ky)=3
             else if (age(kx,ky).eq.5) then
                age(kx,ky)=1
             end if
          end do
         end do  
    
c.     Ingest the forecast from ADV26. In this code, we are ingesting the first "florida-relevent" forecast from Hurricane Irma. 
c.     This file has information on wind, surge, rain risk, and the estimated time of arrival of tropical storm force winds for each grid cell
c.     defined as mx, my. 

       open(unit=26,file='adv26.csv',status='old')
       do while(.true.)
         read(26,*,end=26) mx,my,wrisk(mx,my),srisk(mx,my),
     *        rrisk(mx,my),tstoa(mx,my)
       end do 
 26    close(26)
     
c ... initalize ABM with evacuation-decisions 
c     this is a subroutine at the bottom of the .f file 
       call ABM(1,npop2,wrisk,srisk,rrisk,mevac,pevac,
     *     socio,age,nocar,mobl,agewt,
     *     wwt,swt,rwt,barr,mobwt,strmwt,ewt,tstoa)       
	 
c ... here is where we update the forecast every 6 hours 
       ncount=0
       nmin=0
       ntick=0
       ncella=0
       ntime=0
       ncellb=0
       ntimeb=0
	   nlast=0
       state=0
       nevac=0
	   
c     we loop through every ntick (1.2 seconds).. 
c     forecasts update every 6 hours 
c     there are 3000 nticks in an hour, so every 6 hours is 18000 nticks 
      
        do while(.true.)
         ntick=ntick+1

c.       Ingest the forecast from ADV27     
         IF (ntick.eq.18000) then 
         open(unit=27,file='adv27.csv',status='old')
         do while(.true.)
         read(27,*,end=27) mx,my,wrisk(mx,my),srisk(mx,my),
     *        rrisk(mx,my),tstoa(mx,my)
         end do 
 27     close(27)
         end if   
c.       Ingest the forecast from ADV28    
         IF (ntick.eq.(18000*2)) then 
         open(unit=28,file='adv28.csv',status='old')
         do while(.true.)
         read(28,*,end=28) mx,my,wrisk(mx,my),srisk(mx,my),
     *        rrisk(mx,my),tstoa(mx,my)
         end do 
 28     close(28)
         end if         
c.       Ingest the forecast from ADV29      
         IF (ntick.eq.(18000*3)) then 
         open(unit=29,file='adv29.csv',status='old')
         do while(.true.)
         read(29,*,end=29) mx,my,wrisk(mx,my),srisk(mx,my),
     *        rrisk(mx,my),tstoa(mx,my)
         end do 
 29     close(29)
         end if 
c.       Ingest the forecast from ADV30  
         IF (ntick.eq.(18000*4)) then 
         open(unit=30,file='adv30.csv',status='old')
         do while(.true.)
         read(30,*,end=30) mx,my,wrisk(mx,my),srisk(mx,my),
     *        rrisk(mx,my),tstoa(mx,my)
         end do 
 30     close(30)
         end if 
c.       Ingest the forecast from ADV31  
         IF (ntick.eq.(18000*5)) then 
         open(unit=31,file='adv31.csv',status='old')
         do while(.true.)
         read(31,*,end=31) mx,my,wrisk(mx,my),srisk(mx,my),
     *        rrisk(mx,my),tstoa(mx,my)
         end do 
 31     close(31)
         end if 
c.       Ingest the forecast from ADV32  
         IF (ntick.eq.(18000*6)) then 
         open(unit=32,file='adv32.csv',status='old')
         do while(.true.)
         read(32,*,end=32) mx,my,wrisk(mx,my),srisk(mx,my),
     *        rrisk(mx,my),tstoa(mx,my)
         end do 
 32     close(32)
         end if          
c.       Ingest the forecast from ADV33        
         IF (ntick.eq.(18000*7)) then 
         open(unit=33,file='adv33.csv',status='old')
         do while(.true.)
         read(33,*,end=33) mx,my,wrisk(mx,my),srisk(mx,my),
     *        rrisk(mx,my),tstoa(mx,my)
         end do 
 33     close(33)
         end if
c.       Ingest the forecast from ADV34  
         IF (ntick.eq.(18000*8)) then 
         open(unit=34,file='adv34.csv',status='old')
         do while(.true.)
         read(34,*,end=34) mx,my,wrisk(mx,my),srisk(mx,my),
     *        rrisk(mx,my),tstoa(mx,my)
         end do 
 34     close(34)
         end if  
c.       Ingest the forecast from ADV35  
         IF (ntick.eq.(18000*9)) then 
         open(unit=35,file='adv35.csv',status='old')
         do while(.true.)
         read(35,*,end=35) mx,my,wrisk(mx,my),srisk(mx,my),
     *        rrisk(mx,my),tstoa(mx,my)
         end do 
 35     close(35)
         end if 
c.       Ingest the forecast from ADV36
         IF (ntick.eq.(18000*10)) then 
         open(unit=36,file='adv36.csv',status='old')
         do while(.true.)
         read(36,*,end=36) mx,my,wrisk(mx,my),srisk(mx,my),
     *        rrisk(mx,my),tstoa(mx,my)
         end do 
 36     close(36)
         end if          
c.       Ingest the forecast from ADV37 (appx 0.5 days pre-landfal)       
         IF (ntick.eq.(18000*11)) then 
         open(unit=37,file='adv37.csv',status='old')
         do while(.true.)
         read(37,*,end=37) mx,my,wrisk(mx,my),srisk(mx,my),
     *        rrisk(mx,my),tstoa(mx,my)
         end do 
 37     close(37)
         end if
c.       Ingest the forecast from ADV38      
         IF (ntick.eq.(18000*12)) then 
         open(unit=38,file='adv38.csv',status='old')
         do while(.true.)
         read(38,*,end=38) mx,my,wrisk(mx,my),srisk(mx,my),
     *        rrisk(mx,my),tstoa(mx,my)
         end do 
 38     close(38)
         end if   
c.       Ingest the forecast from ADV39       
         IF (ntick.eq.(18000*13)) then 
         open(unit=39,file='adv39.csv',status='old')
         do while(.true.)
         read(39,*,end=39) mx,my,wrisk(mx,my),srisk(mx,my),
     *        rrisk(mx,my),tstoa(mx,my)
         end do 
 39     close(39)
         end if 
c.       Ingest the forecast from ADV40       
         IF (ntick.eq.(18000*14)) then 
         open(unit=40,file='adv40.csv',status='old')
         do while(.true.)
         read(40,*,end=40) mx,my,wrisk(mx,my),srisk(mx,my),
     *        rrisk(mx,my),tstoa(mx,my)
         end do 
 40     close(40)
         end if         
c.       Ingest the forecast from ADV41 (loc in nw fl)       
         IF (ntick.eq.(18000*15)) then 
         open(unit=41,file='adv41.csv',status='old')
         do while(.true.)
         read(41,*,end=41) mx,my,wrisk(mx,my),srisk(mx,my),
     *        rrisk(mx,my),tstoa(mx,my)
         end do 
 41     close(41)
         end if
c.       Ingest the forecast from ADV42  
         IF (ntick.eq.(18000*16)) then 
         open(unit=42,file='adv42.csv',status='old')
         do while(.true.)
         read(42,*,end=42) mx,my,wrisk(mx,my),srisk(mx,my),
     *        rrisk(mx,my),tstoa(mx,my)
         end do 
 42     close(42)
         end if  
c.       Ingest the forecast from ADV43  
         IF (ntick.eq.(18000*17)) then 
         open(unit=43,file='adv43.csv',status='old')
         do while(.true.)
         read(43,*,end=43) mx,my,wrisk(mx,my),srisk(mx,my),
     *        rrisk(mx,my),tstoa(mx,my)
         end do 
 43     close(43)
         end if 
c.       Ingest the forecast from ADV44 
         IF (ntick.eq.(18000*18)) then 
         open(unit=44,file='adv44.csv',status='old')
         do while(.true.)
         read(44,*,end=44) mx,my,wrisk(mx,my),srisk(mx,my),
     *        rrisk(mx,my),tstoa(mx,my)
         end do 
 44     close(44)
         end if          
c.       Ingest the forecast from ADV45 (appx 0.5 days pre-landfal)       
         IF (ntick.eq.(18000*19)) then 
         open(unit=45,file='adv45.csv',status='old')
         do while(.true.)
         read(45,*,end=45) mx,my,wrisk(mx,my),srisk(mx,my),
     *        rrisk(mx,my),tstoa(mx,my)
         end do 
 45     close(45)
         end if
c.       Ingest the forecast from ADV46       
         IF (ntick.eq.(18000*20)) then 
         open(unit=46,file='adv46.csv',status='old')
         do while(.true.)
         read(46,*,end=46) mx,my,wrisk(mx,my),srisk(mx,my),
     *        rrisk(mx,my),tstoa(mx,my)
         end do 
 46     close(46)
         end if   
c.       Ingest the forecast from ADV47       
         IF (ntick.eq.(18000*21)) then 
         open(unit=47,file='adv47.csv',status='old')
         do while(.true.)
         read(47,*,end=47) mx,my,wrisk(mx,my),srisk(mx,my),
     *        rrisk(mx,my),tstoa(mx,my)
         end do 
 47     close(47)
         end if 
c.       Ingest the forecast from ADV48       
         IF (ntick.eq.(18000*22)) then 
         open(unit=48,file='adv48.csv',status='old')
         do while(.true.)
         read(48,*,end=48) mx,my,wrisk(mx,my),srisk(mx,my),
     *        rrisk(mx,my),tstoa(mx,my)
         end do 
 48     close(48)
         end if 		 

c.  within the same do loop, we have our emergency manager agents make evacuation decisions
c.  this is a sub loop here going across all grid cells 
c.  defining the time when evacuation orders will be issued based on the forecast and clearance times
c.  higher surge risk means and earlier decision point as more time is needed to evacuate  
		 			 		 
         do ky=1,20
         do kx=1,8
          if (ntick.eq.1) then 
            if (srisk(kx,ky).ge.5.and.srisk(kx,ky).le.6) then
              eotime(kx,ky)=int(3000*(tstoa(kx,ky)-ctimey(kx,ky)))
            else if (srisk(kx,ky).ge.7.and.srisk(kx,ky).le.8) then
              eotime(kx,ky)=int(3000*(tstoa(kx,ky)-ctimeo(kx,ky)))
            else if (srisk(kx,ky).ge.9.and.srisk(kx,ky).le.10) then
              eotime(kx,ky)=int(3000*(tstoa(kx,ky)-ctimer(kx,ky)))		  
            else
              eotime(kx,ky)=9999999
            end if   
            
          else if (mod(ntick,18000).eq.0) then
            if (srisk(kx,ky).ge.5.and.srisk(kx,ky).le.6) then
              eotime(kx,ky)=int(3000*(tstoa(kx,ky)-ctimey(kx,ky)))
            else if (srisk(kx,ky).ge.7.and.srisk(kx,ky).le.8) then
              eotime(kx,ky)=int(3000*(tstoa(kx,ky)-ctimeo(kx,ky)))
            else if (srisk(kx,ky).ge.9.and.srisk(kx,ky).le.10) then
              eotime(kx,ky)=int(3000*(tstoa(kx,ky)-ctimer(kx,ky)))            
            else
              eotime(kx,ky)=9999999
            end if 
            if (eotime(kx,ky).lt.0) then
                eotime(kx,ky)=0
            end if                
          end if
         end do
         end do

c.  looping across grid cells again, an evacuation order will be issued once the time (ntick) reaches the evacuation order time (eotime)
         do ky=1,20
         do kx=1,8
         
         if (evac(kx,ky).eq.0) then
            if(ntick.ge.eotime(kx,ky)) then
               evac(kx,ky)=1
            else
               evac(kx,ky)=0
            end if

c.   evacuation orders aren't issued for the entire grid cell, just portions along the coast
c..  this code issues evacuation orders to that percentage of people in the grid cell within an evacuation zone, as defined in the npop file  
      
            IF (evac(kx,ky).eq.1) then
                n3e=npop2(kx,ky)/4 
                do n4=1,n3e  
                   if (rand().le.(parteos(kx,ky)*0.01)) then
                      pevac(kx,ky,n4)=1
                   else 
                      pevac(kx,ky,n4)=0
                   END IF
                end do
            end if 
 
         end if     

         end do
         end do
c.  printing out the evacuation order status every 6 hours
c. this prints to the screen and to different files e.g., fort.10, fort.16         
         
         if (ntick.eq.1) then
              write(*,*) 'evaco time, evac issued'          
              write(16,*) 'evaco time, evac issued'          
           do my=20,1,-1
              write(*,*) my,(eotime(mx,my)/3000,mx=1,8),
     *                     (evac(mx,my),mx=1,8)
              write(10,*) my,(eotime(mx,my)/3000,mx=1,8),
     *                     (evac(mx,my),mx=1,8)     
              write(16,*) my,(eotime(mx,my)/3000,mx=1,8),
     *                     (evac(mx,my),mx=1,8)    
           end do
         end if
         if (mod(ntick,18000).eq.0) then
              write(*,*) 'evaco time, evac issued'  
              write(16,*) 'evaco time, evac issued'  
              write(16,*) 'nhour',ntick/3000    
           do my=20,1,-1
              write(*,*) my,((eotime(mx,my)-ntick)/3000,mx=1,8),
     *                     (evac(mx,my),mx=1,8)
              write(10,*) my,((eotime(mx,my)-ntick)/3000,mx=1,8),
     *                     (evac(mx,my),mx=1,8)
              write(16,*) my,((eotime(mx,my)-ntick)/3000,mx=1,8),
     *                     (evac(mx,my),mx=1,8)    
           end do
         end if         
		 
c.    THE TRAFFIC MODEL
c ... check for any accident/out of gas incident anywhere at this time step
c     incidents have the effect of blocking traffic through the area for 5 minutes
c     and slowing it down for another 10 minutes thereafter

       DO ln=1,5  
c. this section will place random accidents on the E Interstate (if any)
         if (ntick.le.60*50) goto 780
         if (rand().le.rA.and.ntime.eq.0.and.ncella.eq.0) then
           ntrial=0
 777       ncella=1+int(real(NCELLY-1)*rand())
           ntrial=ntrial+1
           if (car(ncella,ln).eq.0.and.ntrial.lt.10) goto 777
           if (car(ncella,ln).eq.0) then
             ncella=0
             goto 778
           end if
           write(*,*) 'E accident@ (ticks,cell,lane): ',ntick,ncella,ln
           write(10,*) 'E accident@ (ticks,cell,lane): ',ntick,ncella,ln
         end if
 778     continue
c. this section will place random accidents on the W Interstate (if any)
         if (rand().le.rA.and.ntimeb.eq.0.and.ncellb.eq.0) then
           ntrial=0
 779       ncellb=1+int(real(NCELLY-1)*rand())
           ntrial=ntrial+1
           if (car(ncellb,ln).eq.0.and.ntrial.lt.10) goto 779
           if (car(ncellb,ln).eq.0) then
             ncellb=0
             goto 780
           end if
           write(*,*) 'W accident@ (ticks,cell,lane): ',ntick,ncellb,ln
           write(10,*) 'W accident@ (ticks,cell,lane): ',ntick,ncellb,ln
         end if
 780     continue
       END DO 
             
c .... MOVING TRAFFIC on the E-Interstate ............
       do ln=1,5
         do i=NCELLY-1,1,-1
c ... Rule0 - accident:
c          incidents have the effect of blocking traffic through the area for 5 minutes 
c          (ntime=500) and slowing it down for another 10 minutes thereafter (ntime=752)
           if (i.eq.ncella.and.ntime.eq.0) ntime=752
           if (ntime.eq.1) then
             ntime=0
             ncella=0
           end if
           if (ntime.gt.0.and.i.eq.ncella) then
             ntime=ntime-1
             nt=0
             if (ntime.le.500.and.car(i+1,ln).eq.0) nt=1
             goto 100
           end if
c ... Rule1 - acceleration:
c       if V<Vmx accel by 1 unless car is stopped
c       if car is stopped, then wait 1 tic before accelerating
           nt=0
           if (car(i,ln).eq.1.and.vel(i,ln).lt.nVMX) then
             if (vel(i,ln).ne.0) then
               nt=min(nVMX,vel(i,ln)+1)
             else
               if (lag(i,ln).eq.1) then
                 lag(i,ln)=0
                 nt=1
               else
                 lag(i,ln)=1
                 nt=0
               end if
             end if
           else if (car(i,ln).eq.1.and.vel(i,ln).eq.nVMX) then
             nt=nVMX
           else
             car2(i,ln)=0
             carp2(i,ln)=0
             vel2(i,ln)=0
             ident2(i,ln)=0
             exit2(i,ln)=0
             two2e(i,ln)=0
             goto 101
           end if
c ... handle if car is at intersection for middle highway
           if (i+nt.gt.17500.and.twoe(i,ln).eq.2) then
             do le=1,3           
               if (carmw(8,10,i+nt-17500,le).eq.0) then
                 carmw(8,10,i+nt-17500,le)=1
                 velmw(8,10,i+nt-17500,le)=vel(i,ln)
                 identmw(8,10,i+nt-17500,le)=ident(i,ln)
                 carpmw(8,10,i+nt-17500,le)=carp(i,ln)
                 twomw(8,10,i+nt-17500,le)=twoe(i,ln)
                 exitmw(8,10,i+nt-17500,le)=exite(i,ln)
                 car2(i,ln)=0
                 carp2(i,ln)=0
                 vel2(i,ln)=0
                 ident2(i,ln)=0
                 exit2(i,ln)=0
                 two2e(i,ln)=0
                 goto 101
               end if              
             end do
                 car2(i,ln)=1
                 vel2(i,ln)=0
                 ident2(i,ln)=ident(i,ln)
                 carp2(i,ln)=carp(i,ln)
                 exit2(i,ln)=exite(i,ln)
                 two2e(i,ln)=twoe(i,ln)
                 goto 101
           end if
c ... handle if car is at edge
           NCELLe=NCELLY
           if (ivac.eq.1) NCELLe=1750*(exit2(i,ln)-1)
           if (i+nt.gt.NCELLe) then
             car2(i,ln)=0
             vel2(i,ln)=0
             ident2(i,ln)=0
             nevac=nevac+carp(i,ln)
             carp2(i,ln)=0
             exit2(i,ln)=0
             two2e(i,ln)=0
             goto 101
           end if
c ... Rule2 - decceleration
c       search from one cell ahead to nt ahead if anywhere in that range there exists a
c       car ahead, then slow down
           i21=i+1
           i22=i+nt
           if (nt.eq.0) goto 100
           do i2=i21,i22
             if (car(i2,ln).eq.1) then
               nt=min(max(0,i2-i-1),3)
               goto 100
             end if
           end do
 100       continue
c ... Rule3 - stochastic
c       randomly slow by one if moving
           if (nt.gt.0.and.rand().le.RP) nt=nt-1
c ... Rule 4 - move cars
c       Move car at its speed to new cell location
           if (car(i,ln).eq.1) then
             if (i+nt.le.NCELLY) then            
               car2(i+nt,ln)=1
               vel2(i+nt,ln)=nt
               ident2(i+nt,ln)=ident(i,ln)
               carp2(i+nt,ln)=carp(i,ln)
               exit2(i+nt,ln)=exite(i,ln)
               two2e(i+nt,ln)=twoe(i,ln)
             end if
             if (nt.gt.0) then
               car2(i,ln)=0
               carp2(i,ln)=0
               vel2(i,ln)=0
               ident2(i,ln)=0
               exit2(i,ln)=0
               two2e(i,ln)=0
             end if
           end if
 101       continue
         end do
       end do   
c .... MOVING TRAFFIC on the W Interstate ............
       do ln=1,5
         do i=NCELLY-1,1,-1
c ... Rule0 - accident:
c          incidents have the effect of blocking traffic through the area for 5 minutes 
c          (ntime=752) and slowing it down for another 10 minutes thereafter (ntime=500)
           if (i.eq.ncellb.and.ntimeb.eq.0) ntimeb=752
           if (ntimeb.eq.1) then
             ntimeb=0
             ncellb=0
           end if
           if (ntimeb.gt.0.and.i.eq.ncellb) then
             ntimeb=ntimeb-1
             nt=0
             if (ntimeb.le.500.and.carX(i+1,ln).eq.0) nt=1
             goto 800
           end if
c ... Rule1 - acceleration:
c       if V<Vmx accel by 1 unless car is stopped
c       if car is stopped, then wait 1 tic before accelerating
           nt=0
           if (carX(i,ln).eq.1.and.velX(i,ln).lt.nVMX) then
             if (velX(i,ln).ne.0) then
               nt=min(nVMX,velX(i,ln)+1)
             else
               if (lagX(i,ln).eq.1) then
                 lagX(i,ln)=0
                 nt=1
               else
                 lagX(i,ln)=1
                 nt=0
               end if
             end if
           else if (carX(i,ln).eq.1.and.velX(i,ln).eq.nVMX) then
             nt=nVMX
           else
             car2X(i,ln)=0
             carp2X(i,ln)=0
             vel2X(i,ln)=0
             ident2X(i,ln)=0
             one2w(i,ln)=0
             two2w(i,ln)=0
             exit2w(i,ln)=0
             goto 801
           end if
c ... handle if car is at intersection for middle highway
           if (i+nt.gt.17500.and.twow(i,ln).eq.1) then
             do le=1,3           
               if (carme(1,10,i+nt-17500,le).eq.0) then
                 carme(1,10,i+nt-17500,le)=1
                 velme(1,10,i+nt-17500,le)=velX(i,ln)
                 identme(1,10,i+nt-17500,le)=identX(i,ln)
                 carpme(1,10,i+nt-17500,le)=carpX(i,ln)
                 twome(1,10,i+nt-17500,le)=twow(i,ln)
                 exitme(1,10,i+nt-17500,le)=exitw(i,ln)
                 car2X(i,ln)=0
                 carp2X(i,ln)=0
                 vel2X(i,ln)=0
                 ident2X(i,ln)=0
                 exit2w(i,ln)=0
                 two2w(i,ln)=0
                 goto 801
               end if              
             end do
                 car2X(i,ln)=1
                 vel2X(i,ln)=0
                 ident2X(i,ln)=identX(i,ln)
                 carp2X(i,ln)=carpX(i,ln)
                 exit2w(i,ln)=exitw(i,ln)
                 two2w(i,ln)=twow(i,ln)
                 goto 801
           end if           
c ... handle if car is at edge
           NCELLe=NCELLY
           if (ivac.eq.1) NCELLe=1750*(exit2w(i,ln)-1)
           if (i+nt.gt.NCELLe) then
             car2X(i,ln)=0
             vel2X(i,ln)=0
             ident2X(i,ln)=0
             nevac=nevac+carpX(i,ln)
             carp2X(i,ln)=0
             one2w(i,ln)=0
             two2w(i,ln)=0
             exit2w(i,ln)=0
             goto 801
           end if
c ... Rule2 - decceleration
c       search from one cell ahead to nt ahead if anywhere in that range there exists a
c       car ahead, then slow down
           i21=i+1
           i22=i+nt
           if (nt.eq.0) goto 800
           do i2=i21,i22
             if (carX(i2,ln).eq.1) then
               nt=min(max(0,i2-i-1),3)
               goto 800
             end if
           end do
 800       continue
c ... Rule3 - stochastic
c       randomly slow by one if moving
           if (nt.gt.0.and.rand().le.RP) nt=nt-1
c ... Rule 4 - move cars
c       Move car at its speed to new cell location
           if (carX(i,ln).eq.1) then
             if (i+nt.le.NCELLY) then            
               car2X(i+nt,ln)=1
               vel2X(i+nt,ln)=nt
               ident2X(i+nt,ln)=identX(i,ln)
               carp2X(i+nt,ln)=carpX(i,ln)
               one2w(i+nt,ln)=onew(i,ln)
               two2w(i+nt,ln)=twow(i,ln)
               exit2w(i+nt,ln)=exitw(i,ln)
             end if
             if (nt.gt.0) then
               car2X(i,ln)=0
               carp2X(i,ln)=0
               vel2X(i,ln)=0
               ident2X(i,ln)=0
               one2w(i,ln)=0
               two2w(i,ln)=0
               exit2w(i,ln)=0
             end if
           end if
 801       continue
         end do
       end do   
       
c .... MOVING TRAFFIC on the S Interstate (eastbound) ............
         ky=1
         do kx=1,8,1
         do le=1,3
         do i=NCELLX,1,-1 
c ... Rule1 - acceleration:
c       if V<Vmx accel by 1 unless car is stopped
c       if car is stopped, then wait 1 tic before accelerating
           nt=0   
           if (carse(kx,ky,i,le).eq.1.and. 
     *           velse(kx,ky,i,le).lt.nVMX) THEN
             if (velse(kx,ky,i,le).ne.0) then
               nt=min(nVMX,velse(kx,ky,i,le)+1)
             else
               if (lagse(kx,ky,i,le).eq.1) then
                 lagse(kx,ky,i,le)=0
                 nt=1
               else
                 lagse(kx,ky,i,le)=1
                 nt=0
               end if
             end if
           else if (carse(kx,ky,i,le).eq.1.and.
     *           velse(kx,ky,i,le).eq.nVMX) then
             nt=nVMX
           else
             car2se(kx,ky,i,le)=0
             carp2se(kx,ky,i,le)=0
             vel2se(kx,ky,i,le)=0
             ident2se(kx,ky,i,le)=0
             one2se(kx,ky,i,le)=0
             two2se(kx,ky,i,le)=0
             exit2se(kx,ky,i,le)=0
             goto 2001
           end if
           if (kx.eq.8) then
c ... handle if car is at edge or moving beyond - move to N/S
             if (i+nt.gt.NCELLX) then
             do ln=1,5           
               if (car2((ky-1)*1750+i-NCELLX+nt,ln).eq.0) then
                 car2((ky-1)*1750+i-NCELLX+nt,ln)=1
                 vel2((ky-1)*1750+i-NCELLX+nt,ln)=velse(kx,ky,i,le)
                 ident2((ky-1)*1750+i-NCELLX+nt,ln)=identse(kx,ky,i,le)
                 carp2((ky-1)*1750+i-NCELLX+nt,ln)=carpse(kx,ky,i,le)
                 rt2((ky-1)*1750+i-NCELLX+nt,ln)=onese(kx,ky,i,le)
                 two2e((ky-1)*1750+i-NCELLX+nt,ln)=twose(kx,ky,i,le)
                 exit2((ky-1)*1750+i-NCELLX+nt,ln)=exitse(kx,ky,i,le)                 
                 car2se(kx,ky,i,le)=0
                 carp2se(kx,ky,i,le)=0
                 vel2se(kx,ky,i,le)=0
                 ident2se(kx,ky,i,le)=0
                 one2se(kx,ky,i,le)=0
                 two2se(kx,ky,i,le)=0
                 exit2se(kx,ky,i,le)=0                  
                 goto 2001
               end if
             end do
                 car2se(kx,ky,i,le)=1
                 vel2se(kx,ky,i,le)=0
                 ident2se(kx,ky,i,le)=identse(kx,ky,i,le)
                 carp2se(kx,ky,i,le)=carpse(kx,ky,i,le)
                 one2se(kx,ky,i,le)=onese(kx,ky,i,le)                   
                 two2se(kx,ky,i,le)=twose(kx,ky,i,le)
                 exit2se(kx,ky,i,le)=exitse(kx,ky,i,le)                 
                 goto 2001
             end if
           else
c ... handle case of car reaching edge of a tile
             if (i+nt.gt.NCELLX) then
               if (car2se(kx+1,ky,i+nt-NCELLX,le).eq.0) then
                 car2se(kx+1,ky,i+nt-NCELLX,le)=1
                 vel2se(kx+1,ky,i+nt-NCELLX,le)=velse(kx,ky,i,le)
                 ident2se(kx+1,ky,i+nt-NCELLX,le)=identse(kx,ky,i,le)
                 carp2se(kx+1,ky,i+nt-NCELLX,le)=carpse(kx,ky,i,le)
                 one2se(kx+1,ky,i+nt-NCELLX,le)=onese(kx,ky,i,le)                         
                 two2se(kx+1,ky,i+nt-NCELLX,le)=twose(kx,ky,i,le)                                 
                 exit2se(kx+1,ky,i+nt-NCELLX,le)=exitse(kx,ky,i,le)                  
                 car2se(kx,ky,i,le)=0
                 carp2se(kx,ky,i,le)=0
                 vel2se(kx,ky,i,le)=0
                 ident2se(kx,ky,i,le)=0
                 one2se(kx,ky,i,le)=0
                 two2se(kx,ky,i,le)=0
                 exit2se(kx,ky,i,le)=0
                 goto 2001
               else
                 car2se(kx,ky,i,le)=1
                 one2se(kx,ky,i,le)=onese(kx,ky,i,le)  
                 vel2se(kx,ky,i,le)=0
                 ident2se(kx,ky,i,le)=identse(kx,ky,i,le)
                 carp2se(kx,ky,i,le)=carpse(kx,ky,i,le)
                 two2se(kx,ky,i,le)=twose(kx,ky,i,le)
                 exit2se(kx,ky,i,le)=exitse(kx,ky,i,le)
                 goto 2001
               end if
             end if
           end if
c ... Rule2 - decceleration
c       search from one cell ahead to nt ahead if anywhere in that range there exists a
c       car ahead, then slow to one cell less to avoid a crash
           i21=i+1
           i22=min(i+nt,NCELLX)
           if (i22.eq.NCELLX) nt=NCELLX-i
           if (nt.eq.0) goto 2000
           do i2=i21,i22
             if (carse(kx,ky,i2,le).eq.1) then
               nt=min(max(0,i2-i-1),1)
               goto 2000
             end if
           end do
 2000      continue
c ... Rule3 - stochastic
c       randomly slow by one if moving
           if (nt.gt.0.and.rand().le.RP) nt=nt-1
c ... Rule 4 - move cars
c       Move car at its speed to new cell location
           if (carse(kx,ky,i,le).eq.1) then
             if (i+nt.le.NCELLX) then            
               car2se(kx,ky,i+nt,le)=1
               one2se(kx,ky,i+nt,le)=onese(kx,ky,i,le)  
               vel2se(kx,ky,i+nt,le)=nt
               ident2se(kx,ky,i+nt,le)=identse(kx,ky,i,le)
               carp2se(kx,ky,i+nt,le)=carpse(kx,ky,i,le)
               two2se(kx,ky,i+nt,le)=twose(kx,ky,i,le)
               exit2se(kx,ky,i+nt,le)=exitse(kx,ky,i,le)
             end if
             if (nt.gt.0) then
               car2se(kx,ky,i,le)=0
               carp2se(kx,ky,i,le)=0
               vel2se(kx,ky,i,le)=0
               ident2se(kx,ky,i,le)=0
               one2se(kx,ky,i,le)=0  
               two2se(kx,ky,i,le)=0
               exit2se(kx,ky,i,le)=0                
             end if
           end if
 2001      continue
         end do
         end do
         end do
         
c .... MOVING TRAFFIC on the S Interstate (westbound) ............
         ky=1
         do kx=8,1,-1
         do le=1,3
         do i=NCELLX,1,-1
c ... Rule1 - acceleration:
c       if V<Vmx accel by 1 unless car is stopped
c       if car is stopped, then wait 1 tic before accelerating
           nt=0
           if (carsw(kx,ky,i,le).eq.1.and. 
     *           velsw(kx,ky,i,le).lt.nVMX) THEN
             if (velsw(kx,ky,i,le).ne.0) then
               nt=min(nVMX,velsw(kx,ky,i,le)+1)
             else
               if (lagsw(kx,ky,i,le).eq.1) then
                 lagsw(kx,ky,i,le)=0
                 nt=1
               else
                 lagsw(kx,ky,i,le)=1
                 nt=0
               end if
             end if
           else if (carsw(kx,ky,i,le).eq.1.and.
     *           velsw(kx,ky,i,le).eq.nVMX) then
             nt=nVMX
           else
             car2sw(kx,ky,i,le)=0
             carp2sw(kx,ky,i,le)=0
             vel2sw(kx,ky,i,le)=0
             ident2sw(kx,ky,i,le)=0
             one2sw(kx,ky,i,le)=0
             two2sw(kx,ky,i,le)=0
             exit2sw(kx,ky,i,le)=0             
             goto 301
           end if
           if (kx.eq.1) then
c ... handle if car is at edge or moving beyond - move to N/S
             if (i+nt.gt.NCELLX) then
             do ln=1,5           
               if (car2X((ky-1)*1750+i-NCELLX+nt,ln).eq.0) then
                 car2X((ky-1)*1750+i-NCELLX+nt,ln)=1
                 vel2X((ky-1)*1750+i-NCELLX+nt,ln)=velsw(kx,ky,i,le)
                 ident2X((ky-1)*1750+i-NCELLX+nt,ln)=identsw(kx,ky,i,le)
                 carp2X((ky-1)*1750+i-NCELLX+nt,ln)=carpsw(kx,ky,i,le)
                 one2w((ky-1)*1750+i-NCELLX+nt,ln)=onesw(kx,ky,i,le)
                 two2w((ky-1)*1750+i-NCELLX+nt,ln)=twosw(kx,ky,i,le)
                 exit2w((ky-1)*1750+i-NCELLX+nt,ln)=exitsw(kx,ky,i,le)
                 car2sw(kx,ky,i,le)=0
                 carp2sw(kx,ky,i,le)=0
                 vel2sw(kx,ky,i,le)=0
                 ident2sw(kx,ky,i,le)=0
                 one2sw(kx,ky,i,le)=0
                 two2sw(kx,ky,i,le)=0
                 exit2sw(kx,ky,i,le)=0 
                 goto 301
               end if
             end do
                 car2sw(kx,ky,i,le)=1
                 vel2sw(kx,ky,i,le)=0
                 ident2sw(kx,ky,i,le)=identsw(kx,ky,i,le)
                 carp2sw(kx,ky,i,le)=carpsw(kx,ky,i,le)
                 one2sw(kx,ky,i,le)=onesw(kx,ky,i,le)
                 two2sw(kx,ky,i,le)=twosw(kx,ky,i,le)
                 exit2sw(kx,ky,i,le)=exitsw(kx,ky,i,le)   
                 goto 301
             end if
           else
c ... handle case of car reaching edge of a tile
             if (i+nt.gt.NCELLX) then
               if (car2sw(kx-1,ky,i+nt-NCELLX,le).eq.0) then
                 car2sw(kx-1,ky,i+nt-NCELLX,le)=1
                 vel2sw(kx-1,ky,i+nt-NCELLX,le)=velsw(kx,ky,i,le)
                 ident2sw(kx-1,ky,i+nt-NCELLX,le)=identsw(kx,ky,i,le)
                 carp2sw(kx-1,ky,i+nt-NCELLX,le)=carpsw(kx,ky,i,le)
                 one2sw(kx-1,ky,i+nt-NCELLX,le)=onesw(kx,ky,i,le)                         
                 two2sw(kx-1,ky,i+nt-NCELLX,le)=twosw(kx,ky,i,le)                                 
                 exit2sw(kx-1,ky,i+nt-NCELLX,le)=exitsw(kx,ky,i,le)                  
                 car2sw(kx,ky,i,le)=0
                 carp2sw(kx,ky,i,le)=0
                 vel2sw(kx,ky,i,le)=0
                 ident2sw(kx,ky,i,le)=0
                 one2sw(kx,ky,i,le)=0
                 two2sw(kx,ky,i,le)=0
                 exit2sw(kx,ky,i,le)=0                 
                 goto 301
               else
                 car2sw(kx,ky,i,le)=1
                 vel2sw(kx,ky,i,le)=0
                 ident2sw(kx,ky,i,le)=identsw(kx,ky,i,le)
                 carp2sw(kx,ky,i,le)=carpsw(kx,ky,i,le)
                 one2sw(kx,ky,i,le)=onesw(kx,ky,i,le)                   
                 two2sw(kx,ky,i,le)=twosw(kx,ky,i,le)
                 exit2sw(kx,ky,i,le)=exitsw(kx,ky,i,le)                 
                 goto 301
               end if
             end if
           end if
c ... Rule2 - decceleration
c       search from one cell ahead to nt ahead if anywhere in that range there exists a
c       car ahead, then slow to one cell less to avoid a crash
           i21=i+1
           i22=min(i+nt,NCELLX)
           if (i22.eq.NCELLX) nt=NCELLX-i
           if (nt.eq.0) goto 300
           do i2=i21,i22
             if (carsw(kx,ky,i2,le).eq.1) then
               nt=min(max(0,i2-i-1),1)
               goto 300
             end if
           end do
 300      continue
c ... Rule3 - stochastic
c       randomly slow by one if moving
           if (nt.gt.0.and.rand().le.RP) nt=nt-1
c ... Rule 4 - move cars
c       Move car at its speed to new cell location
           if (carsw(kx,ky,i,le).eq.1) then
             if (i+nt.le.NCELLX) then            
               car2sw(kx,ky,i+nt,le)=1
               vel2sw(kx,ky,i+nt,le)=nt
               ident2sw(kx,ky,i+nt,le)=identsw(kx,ky,i,le)
               carp2sw(kx,ky,i+nt,le)=carpsw(kx,ky,i,le)
               one2sw(kx,ky,i+nt,le)=onesw(kx,ky,i,le)                 
               two2sw(kx,ky,i+nt,le)=twosw(kx,ky,i,le)
               exit2sw(kx,ky,i+nt,le)=exitsw(kx,ky,i,le)
             end if
             if (nt.gt.0) then
               car2sw(kx,ky,i,le)=0
               carp2sw(kx,ky,i,le)=0
               vel2sw(kx,ky,i,le)=0
               ident2sw(kx,ky,i,le)=0
               one2sw(kx,ky,i,le)=0  
               two2sw(kx,ky,i,le)=0
               exit2sw(kx,ky,i,le)=0                   
             end if
           end if
 301      continue
         end do
         end do
         end do            

c .... MOVING TRAFFIC on the mid-interstate (westbound) ............
         ky=10
         do kx=8,5,-1
         do le=1,3
         do i=NCELLX,1,-1
c ... Rule1 - acceleration:
c       if V<Vmx accel by 1 unless car is stopped
c       if car is stopped, then wait 1 tic before accelerating
           nt=0
           if (carmw(kx,ky,i,le).eq.1.and. 
     *           velmw(kx,ky,i,le).lt.nVMX) THEN
             if (velmw(kx,ky,i,le).ne.0) then
               nt=min(nVMX,velmw(kx,ky,i,le)+1)
             else
               if (lagmw(kx,ky,i,le).eq.1) then
                 lagmw(kx,ky,i,le)=0
                 nt=1
               else
                 lagmw(kx,ky,i,le)=1
                 nt=0
               end if
             end if
           else if (carmw(kx,ky,i,le).eq.1.and.
     *           velmw(kx,ky,i,le).eq.nVMX) then
             nt=nVMX
           else
             car2mw(kx,ky,i,le)=0
             carp2mw(kx,ky,i,le)=0
             vel2mw(kx,ky,i,le)=0
             ident2mw(kx,ky,i,le)=0
             one2mw(kx,ky,i,le)=0
             two2mw(kx,ky,i,le)=0
             exit2mw(kx,ky,i,le)=0             
             goto 3001
           end if
c ... handle if car is at edge
           if (kx.eq.3) then
             if (i+nt.gt.NCELLX) then
             car2mw(kx,ky,i,le)=0
             vel2mw(kx,ky,i,le)=0
             ident2mw(kx,ky,i,le)=0
             nevac=nevac+carpmw(kx,ky,i,le)
             carp2mw(kx,ky,i,le)=0
             exit2mw(kx,ky,i,le)=0   
             two2mw(kx,ky,i,le)=0
             goto 3001
             end if
           end if 
           if (kx.eq.4) then
c ... handle case of car reaching edge of a tile
             if (i+nt.gt.NCELLX) then
               if (carmw(kx-1,ky,i+nt-NCELLX,le).eq.0) then
                 carmw(kx-1,ky,i+nt-NCELLX,le)=1  
                 velmw(kx-1,ky,i+nt-NCELLX,le)=velmw(kx,ky,i,le)
                 identmw(kx-1,ky,i+nt-NCELLX,le)=identmw(kx,ky,i,le)
                 carpmw(kx-1,ky,i+nt-NCELLX,le)=carpmw(kx,ky,i,le)
                 onemw(kx-1,ky,i+nt-NCELLX,le)=onemw(kx,ky,i,le)                         
                 twomw(kx-1,ky,i+nt-NCELLX,le)=twomw(kx,ky,i,le)                                 
                 exitmw(kx-1,ky,i+nt-NCELLX,le)=exitmw(kx,ky,i,le)                          
                 car2mw(kx,ky,i,le)=0
                 carp2mw(kx,ky,i,le)=0
                 vel2mw(kx,ky,i,le)=0
                 ident2mw(kx,ky,i,le)=0
                 one2mw(kx,ky,i,le)=0
                 two2mw(kx,ky,i,le)=0
                 exit2mw(kx,ky,i,le)=0                 
                 goto 3001
               else
                 car2mw(kx,ky,i,le)=1
                 vel2mw(kx,ky,i,le)=0
                 ident2mw(kx,ky,i,le)=identmw(kx,ky,i,le)
                 carp2mw(kx,ky,i,le)=carpmw(kx,ky,i,le)
                 one2mw(kx,ky,i,le)=onemw(kx,ky,i,le)                   
                 two2mw(kx,ky,i,le)=twomw(kx,ky,i,le)
                 exit2mw(kx,ky,i,le)=exitmw(kx,ky,i,le)                 
                 goto 3001
               end if
             end if
           end if
c ... Rule2 - decceleration
c       search from one cell ahead to nt ahead if anywhere in that range there exists a
c       car ahead, then slow to one cell less to avoid a crash
           i21=i+1
           i22=min(i+nt,NCELLX)
           if (i22.eq.NCELLX) nt=NCELLX-i
           if (nt.eq.0) goto 3000
           do i2=i21,i22
             if (carmw(kx,ky,i2,le).eq.1) then
               nt=min(max(0,i2-i-1),1)
               goto 3000
             end if
           end do
 3000      continue
c ... Rule3 - stochastic
c       randomly slow by one if moving
           if (nt.gt.0.and.rand().le.RP) nt=nt-1
c ... Rule 4 - move cars
c       Move car at its speed to new cell location 
           if (carmw(kx,ky,i,le).eq.1) then
             if (i+nt.le.NCELLX) then            
               car2mw(kx,ky,i+nt,le)=1
               vel2mw(kx,ky,i+nt,le)=nt
               ident2mw(kx,ky,i+nt,le)=identmw(kx,ky,i,le)
               carp2mw(kx,ky,i+nt,le)=carpmw(kx,ky,i,le)
               one2mw(kx,ky,i+nt,le)=onemw(kx,ky,i,le)                 
               two2mw(kx,ky,i+nt,le)=twomw(kx,ky,i,le)
               exit2mw(kx,ky,i+nt,le)=exitmw(kx,ky,i,le)
             end if
             if (nt.gt.0) then
               car2mw(kx,ky,i,le)=0
               carp2mw(kx,ky,i,le)=0
               vel2mw(kx,ky,i,le)=0
               ident2mw(kx,ky,i,le)=0
               one2mw(kx,ky,i,le)=0  
               two2mw(kx,ky,i,le)=0
               exit2mw(kx,ky,i,le)=0                   
             end if
           end if
 3001      continue
         end do
         end do
         end do 
c. .... MOVING TRAFFIC on the mid-interstate (eastbound) ............
         ky=10
         do kx=1,4
         do le=1,3
         do i=NCELLX,1,-1
c ... Rule1 - acceleration:
c       if V<Vmx accel by 1 unless car is stopped
c       if car is stopped, then wait 1 tic before accelerating
           nt=0
           if (carme(kx,ky,i,le).eq.1.and. 
     *           velme(kx,ky,i,le).lt.nVMX) THEN
             if (velme(kx,ky,i,le).ne.0) then
               nt=min(nVMX,velme(kx,ky,i,le)+1)
             else
               if (lagme(kx,ky,i,le).eq.1) then
                 lagme(kx,ky,i,le)=0
                 nt=1
               else
                 lagme(kx,ky,i,le)=1
                 nt=0
               end if
             end if
           else if (carme(kx,ky,i,le).eq.1.and.
     *           velme(kx,ky,i,le).eq.nVMX) then
             nt=nVMX
           else
             car2me(kx,ky,i,le)=0
             carp2me(kx,ky,i,le)=0
             vel2me(kx,ky,i,le)=0
             ident2me(kx,ky,i,le)=0
             one2me(kx,ky,i,le)=0
             two2me(kx,ky,i,le)=0
             exit2me(kx,ky,i,le)=0             
             goto 9001
           end if
c ... handle if car is at edge
           if (kx.eq.2) then
             if (i+nt.gt.NCELLX) then
             car2me(kx,ky,i,le)=0
             vel2me(kx,ky,i,le)=0
             ident2me(kx,ky,i,le)=0
             nevac=nevac+carpme(kx,ky,i,le)
             carp2me(kx,ky,i,le)=0
             exit2me(kx,ky,i,le)=0   
             two2me(kx,ky,i,le)=0
             goto 9001
             end if
           end if 
           if (kx.eq.1) then
c ... handle case of car reaching edge of a tile
             if (i+nt.gt.NCELLX) then
               if (carme(kx+1,ky,i+nt-NCELLX,le).eq.0) then
                 carme(kx+1,ky,i+nt-NCELLX,le)=1  
                 velme(kx+1,ky,i+nt-NCELLX,le)=velme(kx,ky,i,le)
                 identme(kx+1,ky,i+nt-NCELLX,le)=identme(kx,ky,i,le)
                 carpme(kx+1,ky,i+nt-NCELLX,le)=carpme(kx,ky,i,le)
                 oneme(kx+1,ky,i+nt-NCELLX,le)=oneme(kx,ky,i,le)                         
                 twome(kx+1,ky,i+nt-NCELLX,le)=twome(kx,ky,i,le)                                 
                 exitme(kx+1,ky,i+nt-NCELLX,le)=exitme(kx,ky,i,le)                          
                 car2me(kx,ky,i,le)=0
                 carp2me(kx,ky,i,le)=0
                 vel2me(kx,ky,i,le)=0
                 ident2me(kx,ky,i,le)=0
                 one2me(kx,ky,i,le)=0
                 two2me(kx,ky,i,le)=0
                 exit2me(kx,ky,i,le)=0                 
                 goto 9001
               else
                 car2me(kx,ky,i,le)=1
                 vel2me(kx,ky,i,le)=0
                 ident2me(kx,ky,i,le)=identme(kx,ky,i,le)
                 carp2me(kx,ky,i,le)=carpme(kx,ky,i,le)
                 one2me(kx,ky,i,le)=oneme(kx,ky,i,le)                   
                 two2me(kx,ky,i,le)=twome(kx,ky,i,le)
                 exit2me(kx,ky,i,le)=exitme(kx,ky,i,le)                 
                 goto 9001
               end if
             end if
           end if
c ... Rule2 - decceleration
c       search from one cell ahead to nt ahead if anywhere in that range there exists a
c       car ahead, then slow to one cell less to avoid a crash
           i21=i+1
           i22=min(i+nt,NCELLX)
           if (i22.eq.NCELLX) nt=NCELLX-i
           if (nt.eq.0) goto 9000
           do i2=i21,i22
             if (carme(kx,ky,i2,le).eq.1) then
               nt=min(max(0,i2-i-1),1)
               goto 9000
             end if
           end do
 9000      continue
c ... Rule3 - stochastic
c       randomly slow by one if moving
           if (nt.gt.0.and.rand().le.RP) nt=nt-1
c ... Rule 4 - move cars
c       Move car at its speed to new cell location 
           if (carme(kx,ky,i,le).eq.1) then
             if (i+nt.le.NCELLX) then            
               car2me(kx,ky,i+nt,le)=1
               vel2me(kx,ky,i+nt,le)=nt
               ident2me(kx,ky,i+nt,le)=identme(kx,ky,i,le)
               carp2me(kx,ky,i+nt,le)=carpme(kx,ky,i,le)
               one2me(kx,ky,i+nt,le)=oneme(kx,ky,i,le)                 
               two2me(kx,ky,i+nt,le)=twome(kx,ky,i,le)
               exit2me(kx,ky,i+nt,le)=exitme(kx,ky,i,le)
             end if
             if (nt.gt.0) then
               car2me(kx,ky,i,le)=0
               carp2me(kx,ky,i,le)=0
               vel2me(kx,ky,i,le)=0
               ident2me(kx,ky,i,le)=0
               one2me(kx,ky,i,le)=0  
               two2me(kx,ky,i,le)=0
               exit2me(kx,ky,i,le)=0                   
             end if
           end if
 9001      continue
         end do
         end do
         end do
c ... MOVING TRAFFIC on the east (2 lane) highways 
         do ky=2,20
		 if (ky.eq.2) goto 6099
		 if (ky.eq.4) goto 6099
		 if (ky.eq.6) goto 6099
		 if (ky.eq.8) goto 6099
		 if (ky.eq.10) goto 6099
		 if (ky.eq.12) goto 6099
		 if (ky.eq.14) goto 6099
		 if (ky.eq.16) goto 6099
		 if (ky.eq.18) goto 6099
		 if (ky.eq.20) goto 6099
         do kx=5,8,1
	     do lx=1,2
         do i=NCELLX,1,-1
         
c ... Rule1 - acceleration:
c       if V<Vmx2 accel by 1 unless car is stopped
c       if car is stopped, then wait 1 tic before accelerating
           nt=0
           if (care(kx,ky,i,lx).eq.1.and.
     *        vele(kx,ky,i,lx).lt.nVMX2) THEN 
             if (vele(kx,ky,i,lx).ne.0) then
               nt=min(nVMX2,vele(kx,ky,i,lx)+1)
             else
               if (lage(kx,ky,i,lx).eq.1) then
                 lage(kx,ky,i,lx)=0
                 nt=1
               else
                 lage(kx,ky,i,lx)=1
                 nt=0
               end if
             end if
           else if (care(kx,ky,i,lx).eq.1.and.
     *        vele(kx,ky,i,lx).eq.nVMX2) then
             nt=nVMX2
           else
             car2e(kx,ky,i,lx)=0
             carp2e(kx,ky,i,lx)=0
             vel2e(kx,ky,i,lx)=0
             ident2e(kx,ky,i,lx)=0
             one2h(kx,ky,i,lx)=0
             two2h(kx,ky,i,lx)=0
             exit2h(kx,ky,i,lx)=0
             goto 1001
           end if
           if (kx.eq.8) then
c ... handle if car is at edge or moving beyond - move to N/S
             if (i+nt.gt.NCELLX) then
             do ln=1,5           
               if (car2((ky-1)*1750+i-NCELLX+nt,ln).eq.0) then
                 car2((ky-1)*1750+i-NCELLX+nt,ln)=1
                 vel2((ky-1)*1750+i-NCELLX+nt,ln)=vele(kx,ky,i,lx)
                 ident2((ky-1)*1750+i-NCELLX+nt,ln)=idente(kx,ky,i,lx)
                 carp2((ky-1)*1750+i-NCELLX+nt,ln)=carpe(kx,ky,i,lx)
                 rt2((ky-1)*1750+i-NCELLX+nt,ln)=oneh(kx,ky,i,lx)
                 two2e((ky-1)*1750+i-NCELLX+nt,ln)=twoh(kx,ky,i,lx)
                 exit2((ky-1)*1750+i-NCELLX+nt,ln)=exith(kx,ky,i,lx) 
                 car2e(kx,ky,i,lx)=0
                 carp2e(kx,ky,i,lx)=0
                 vel2e(kx,ky,i,lx)=0
                 ident2e(kx,ky,i,lx)=0
                 one2h(kx,ky,i,lx)=0
                 two2h(kx,ky,i,lx)=0
                 exit2h(kx,ky,i,lx)=0                  
                 goto 1001
               end if
             end do
                 car2e(kx,ky,i,lx)=1
                 vel2e(kx,ky,i,lx)=0
                 ident2e(kx,ky,i,lx)=idente(kx,ky,i,lx)
                 carp2e(kx,ky,i,lx)=carpe(kx,ky,i,lx)
                 one2h(kx,ky,i,lx)=oneh(kx,ky,i,lx)  
                 two2h(kx,ky,i,lx)=twoh(kx,ky,i,lx)
                 exit2h(kx,ky,i,lx)=exith(kx,ky,i,lx)                  
                 goto 1001
             end if
           else
c ... handle case of car reaching edge of a tile
             if (i+nt.gt.NCELLX) then
               if (car2e(kx+1,ky,i+nt-NCELLX,lx).eq.0) then
                 car2e(kx+1,ky,i+nt-NCELLX,lx)=1
                 vel2e(kx+1,ky,i+nt-NCELLX,lx)=vele(kx,ky,i,lx)
                 ident2e(kx+1,ky,i+nt-NCELLX,lx)=idente(kx,ky,i,lx)
                 carp2e(kx+1,ky,i+nt-NCELLX,lx)=carpe(kx,ky,i,lx)
                 one2h(kx+1,ky,i+nt-NCELLX,lx)=oneh(kx,ky,i,lx)                         
                 two2h(kx+1,ky,i+nt-NCELLX,lx)=twoh(kx,ky,i,lx)                                 
                 exit2h(kx+1,ky,i+nt-NCELLX,lx)=exith(kx,ky,i,lx)                 
                 car2e(kx,ky,i,lx)=0
                 carp2e(kx,ky,i,lx)=0
                 vel2e(kx,ky,i,lx)=0
                 ident2e(kx,ky,i,lx)=0
                 one2h(kx,ky,i,lx)=0
                 two2h(kx,ky,i,lx)=0
                 exit2h(kx,ky,i,lx)=0                 
                 goto 1001
               else
                 car2e(kx,ky,i,lx)=1
                 vel2e(kx,ky,i,lx)=0
                 ident2e(kx,ky,i,lx)=idente(kx,ky,i,lx)
                 carp2e(kx,ky,i,lx)=carpe(kx,ky,i,lx)
                 one2h(kx,ky,i,lx)=oneh(kx,ky,i,lx)  
                 two2h(kx,ky,i,lx)=twoh(kx,ky,i,lx)
                 exit2h(kx,ky,i,lx)=exith(kx,ky,i,lx)
                 goto 1001
               end if
             end if
           end if
c ... Rule2 - decceleration
c       search from one cell ahead to nt ahead if anywhere in that range there exists a
c       car ahead, then slow to one cell less to avoid a crash
           i21=i+1
           i22=min(i+nt,NCELLX)
           if (i22.eq.NCELLX) nt=NCELLX-i
           if (nt.eq.0) goto 1000
           do i2=i21,i22
             if (care(kx,ky,i2,lx).eq.1) then
               nt=min(max(0,i2-i-1),1)
               goto 1000
             end if
           end do
 1000      continue
c ... Rule3 - stochastic
c       randomly slow by one if moving
           if (nt.gt.0.and.rand().le.RP) nt=nt-1
c ... Rule 4 - move cars
c       Move car at its speed to new cell location
           if (care(kx,ky,i,lx).eq.1) then
             if (i+nt.le.NCELLX) then            
               car2e(kx,ky,i+nt,lx)=1
               vel2e(kx,ky,i+nt,lx)=nt
               ident2e(kx,ky,i+nt,lx)=idente(kx,ky,i,lx)
               carp2e(kx,ky,i+nt,lx)=carpe(kx,ky,i,lx)
               one2h(kx,ky,i+nt,lx)=oneh(kx,ky,i,lx)  
               two2h(kx,ky,i+nt,lx)=twoh(kx,ky,i,lx)
               exit2h(kx,ky,i+nt,lx)=exith(kx,ky,i,lx)
             end if
             if (nt.gt.0) then
               car2e(kx,ky,i,lx)=0
               carp2e(kx,ky,i,lx)=0
               vel2e(kx,ky,i,lx)=0
               ident2e(kx,ky,i,lx)=0
               one2h(kx,ky,i,lx)=0  
               two2h(kx,ky,i,lx)=0
               exit2h(kx,ky,i,lx)=0
             end if
           end if
 1001      continue 
         
		 end do
         end do
         end do 
 6099    continue
		 end do

 
c ... MOVING TRAFFIC on the west moving highways 
         do ky=2,20
		 if (ky.eq.2) goto 6098
		 if (ky.eq.4) goto 6098
		 if (ky.eq.6) goto 6098
		 if (ky.eq.8) goto 6098
		 if (ky.eq.10) goto 6098
		 if (ky.eq.12) goto 6098
		 if (ky.eq.14) goto 6098
		 if (ky.eq.16) goto 6098
		 if (ky.eq.18) goto 6098
		 if (ky.eq.20) goto 6098
         do kx=4,1,-1
	     do lx=1,2
         do i=NCELLX,1,-1
c ... Rule1 - acceleration:
c       if V<Vmx2 accel by 1 unless car is stopped
c       if car is stopped, then wait 1 tic before accelerating
           nt=0
           if (care(kx,ky,i,lx).eq.1.and.
     *          vele(kx,ky,i,lx).lt.nVMX2) then
             if (vele(kx,ky,i,lx).ne.0) then
               nt=min(nVMX2,vele(kx,ky,i,lx)+1)
             else
               if (lage(kx,ky,i,lx).eq.1) then
                 lage(kx,ky,i,lx)=0
                 nt=1
               else
                 lage(kx,ky,i,lx)=1
                 nt=0
               end if
             end if
           else if (care(kx,ky,i,lx).eq.1.and.
     *           vele(kx,ky,i,lx).eq.nVMX2) then
             nt=nVMX2
           else
             car2e(kx,ky,i,lx)=0
             carp2e(kx,ky,i,lx)=0
             vel2e(kx,ky,i,lx)=0
             ident2e(kx,ky,i,lx)=0
             one2h(kx,ky,i,lx)=0
             two2h(kx,ky,i,lx)=0
             exit2h(kx,ky,i,lx)=0
             goto 6001
           end if
           if (kx.eq.1) then
c ... handle if car is at edge or moving beyond - move to N/S
             if (i+nt.gt.NCELLX) then
             do ln=1,5           
               if (car2X((ky-1)*1750+i-NCELLX+nt,ln).eq.0) then
                 car2X((ky-1)*1750+i-NCELLX+nt,ln)=1
                 vel2X((ky-1)*1750+i-NCELLX+nt,ln)=vele(kx,ky,i,lx)
                 ident2X((ky-1)*1750+i-NCELLX+nt,ln)=idente(kx,ky,i,lx)
                 carp2X((ky-1)*1750+i-NCELLX+nt,ln)=carpe(kx,ky,i,lx)
                 one2w((ky-1)*1750+i-NCELLX+nt,ln)=oneh(kx,ky,i,lx)
                 two2w((ky-1)*1750+i-NCELLX+nt,ln)=twoh(kx,ky,i,lx)
                 exit2w((ky-1)*1750+i-NCELLX+nt,ln)=exith(kx,ky,i,lx)
                 car2e(kx,ky,i,lx)=0
                 carp2e(kx,ky,i,lx)=0
                 vel2e(kx,ky,i,lx)=0
                 ident2e(kx,ky,i,lx)=0
                 one2h(kx,ky,i,lx)=0
                 two2h(kx,ky,i,lx)=0
                 exit2h(kx,ky,i,lx)=0                 
                 goto 6001
               end if
             end do
              car2e(kx,ky,i,lx)=1
              vel2e(kx,ky,i,lx)=0
              ident2e(kx,ky,i,lx)=idente(kx,ky,i,lx)
              carp2e(kx,ky,i,lx)=carpe(kx,ky,i,lx)
              one2h(kx,ky,i,lx)=oneh(kx,ky,i,lx)  
              two2h(kx,ky,i,lx)=twoh(kx,ky,i,lx)
              exit2h(kx,ky,i,lx)=exith(kx,ky,i,lx)                
              goto 6001 
             end if
           else
c ... handle case of car reaching edge of a tile
             if (i+nt.gt.NCELLX) then
               if (care(kx-1,ky,i+nt-NCELLX,lx).eq.0) then
                 car2e(kx-1,ky,i+nt-NCELLX,lx)=1
                 vel2e(kx-1,ky,i+nt-NCELLX,lx)=vele(kx,ky,i,lx)
                 ident2e(kx-1,ky,i+nt-NCELLX,lx)=idente(kx,ky,i,lx)
                 carp2e(kx-1,ky,i+nt-NCELLX,lx)=carpe(kx,ky,i,lx)
                 one2h(kx-1,ky,i+nt-NCELLX,lx)=oneh(kx,ky,i,lx)                         
                 two2h(kx-1,ky,i+nt-NCELLX,lx)=twoh(kx,ky,i,lx)                                 
                 exit2h(kx-1,ky,i+nt-NCELLX,lx)=exith(kx,ky,i,lx)
                 car2e(kx,ky,i,lx)=0
                 carp2e(kx,ky,i,lx)=0
                 vel2e(kx,ky,i,lx)=0
                 ident2e(kx,ky,i,lx)=0
                 one2h(kx,ky,i,lx)=0
                 two2h(kx,ky,i,lx)=0
                 exit2h(kx,ky,i,lx)=0 
              goto 6001
               else
                 car2e(kx,ky,i,lx)=1
                 vel2e(kx,ky,i,lx)=0
                 ident2e(kx,ky,i,lx)=idente(kx,ky,i,lx)
                 carp2e(kx,ky,i,lx)=carpe(kx,ky,i,lx)
                 one2h(kx,ky,i,lx)=oneh(kx,ky,i,lx)  
                 two2h(kx,ky,i,lx)=twoh(kx,ky,i,lx)
                 exit2h(kx,ky,i,lx)=exith(kx,ky,i,lx)                 
              goto 6001
               end if
             end if
           end if
c ... Rule2 - decceleration
c       search from one cell ahead to nt ahead if anywhere in that range there exists a
c       car ahead, then slow to one cell less to avoid a crash
           i21=i+1
           i22=min(i+nt,NCELLX)
           if (i22.eq.NCELLX) nt=NCELLX-i
           if (nt.eq.0) goto 6000
           do i2=i21,i22
             if (care(kx,ky,i2,lx).eq.1) then
               nt=min(max(0,i2-i-1),1)
               goto 6000
             end if
           end do
 6000      continue
c ... Rule3 - stochastic
c       randomly slow by one if moving
           if (nt.gt.0.and.rand().le.RP) nt=nt-1
c ... Rule 4 - move cars
c       Move car at its speed to new cell location
           if (care(kx,ky,i,lx).eq.1) then
             if (i+nt.le.NCELLX) then            
               car2e(kx,ky,i+nt,lx)=1
               vel2e(kx,ky,i+nt,lx)=nt
               ident2e(kx,ky,i+nt,lx)=idente(kx,ky,i,lx)
               carp2e(kx,ky,i+nt,lx)=carpe(kx,ky,i,lx)
               one2h(kx,ky,i+nt,lx)=oneh(kx,ky,i,lx)  
               two2h(kx,ky,i+nt,lx)=twoh(kx,ky,i,lx)
               exit2h(kx,ky,i+nt,lx)=exith(kx,ky,i,lx)
             end if
             if (nt.gt.0) then
               car2e(kx,ky,i,lx)=0
               carp2e(kx,ky,i,lx)=0
               vel2e(kx,ky,i,lx)=0
               ident2e(kx,ky,i,lx)=0
               one2h(kx,ky,i,lx)=0  
               two2h(kx,ky,i,lx)=0
               exit2h(kx,ky,i,lx)=0
             end if
           end if
 6001      continue
         end do
         end do
         end do
 6098    continue		 
		 end do
c ... now update all cells
         do ln=1,5
         do i=1,NCELLY
           car(i,ln)=car2(i,ln)
           carp(i,ln)=carp2(i,ln)
           vel(i,ln)=vel2(i,ln)
           ident(i,ln)=ident2(i,ln)
           exite(i,ln)=exit2(i,ln)
           twoe(i,ln)=two2e(i,ln)
           carX(i,ln)=car2X(i,ln)
           carpX(i,ln)=carp2X(i,ln)
           velX(i,ln)=vel2X(i,ln)
           identX(i,ln)=ident2X(i,ln)
           onew(i,ln)=one2w(i,ln)
           twow(i,ln)=two2w(i,ln)
           exitw(i,ln)=exit2w(i,ln)
         end do
         end do    
         ky=1
         do kx=1,NX
         do le=1,3
           do i=1,NCELLX  
             carse(kx,ky,i,le)=car2se(kx,ky,i,le)
             carpse(kx,ky,i,le)=carp2se(kx,ky,i,le)
             velse(kx,ky,i,le)=vel2se(kx,ky,i,le)
             identse(kx,ky,i,le)=ident2se(kx,ky,i,le)
             onese(kx,ky,i,le)=one2se(kx,ky,i,le) 
             twose(kx,ky,i,le)=two2se(kx,ky,i,le)
             exitse(kx,ky,i,le)=exit2se(kx,ky,i,le)                   
           end do
         end do
         end DO  
         ky=1
         do kx=8,1,-1
         do le=1,3
           do i=1,NCELLX
             carsw(kx,ky,i,le)=car2sw(kx,ky,i,le)
             carpsw(kx,ky,i,le)=carp2sw(kx,ky,i,le)
             velsw(kx,ky,i,le)=vel2sw(kx,ky,i,le)
             identsw(kx,ky,i,le)=ident2sw(kx,ky,i,le)
             onesw(kx,ky,i,le)=one2sw(kx,ky,i,le) 
             twosw(kx,ky,i,le)=two2sw(kx,ky,i,le)
             exitsw(kx,ky,i,le)=exit2sw(kx,ky,i,le)              
           end do
         end do
         end DO
         ky=10
         do kx=8,5,-1
         do le=1,3
           do i=1,NCELLX
             carmw(kx,ky,i,le)=car2mw(kx,ky,i,le)
             carpmw(kx,ky,i,le)=carp2mw(kx,ky,i,le)
             velmw(kx,ky,i,le)=vel2mw(kx,ky,i,le)
             identmw(kx,ky,i,le)=ident2mw(kx,ky,i,le)
             onemw(kx,ky,i,le)=one2mw(kx,ky,i,le) 
             twomw(kx,ky,i,le)=two2mw(kx,ky,i,le)
             exitmw(kx,ky,i,le)=exit2mw(kx,ky,i,le)              
           end do
         end do
         end DO
         ky=10
         do kx=1,4
         do le=1,3
           do i=1,NCELLX
             carme(kx,ky,i,le)=car2me(kx,ky,i,le)
             carpme(kx,ky,i,le)=carp2me(kx,ky,i,le)
             velme(kx,ky,i,le)=vel2me(kx,ky,i,le)
             identme(kx,ky,i,le)=ident2me(kx,ky,i,le)
             oneme(kx,ky,i,le)=one2me(kx,ky,i,le) 
             twome(kx,ky,i,le)=two2me(kx,ky,i,le)
             exitme(kx,ky,i,le)=exit2me(kx,ky,i,le) 
           end do
         end do
         end DO 
         do lx=1,2        
         do ky=2,NY
		 if (ky.eq.2) goto 6097
		 if (ky.eq.4) goto 6097
		 if (ky.eq.6) goto 6097
		 if (ky.eq.8) goto 6097
		 if (ky.eq.10) goto 6097
		 if (ky.eq.12) goto 6097
		 if (ky.eq.14) goto 6097
		 if (ky.eq.16) goto 6097
		 if (ky.eq.18) goto 6097
		 if (ky.eq.20) goto 6097
         do kx=1,NX
           do i=1,NCELLX   
             care(kx,ky,i,lx)=car2e(kx,ky,i,lx)
             carpe(kx,ky,i,lx)=carp2e(kx,ky,i,lx)
             vele(kx,ky,i,lx)=vel2e(kx,ky,i,lx)
             idente(kx,ky,i,lx)=ident2e(kx,ky,i,lx)
             oneh(kx,ky,i,lx)=one2h(kx,ky,i,lx) 
             twoh(kx,ky,i,lx)=two2h(kx,ky,i,lx)
             exith(kx,ky,i,lx)=exit2h(kx,ky,i,lx)  
           end do
         end do
 6097    continue		 	 
         end do
		 end do
		 
c ... END OF TRAFFIC MODEL SECTION ..............................................c
         
c. Here we are defining some initial parameter states for the evacuation decisions and traffic models		 
c. Jinit is referred to in the subroutine, where it is used to set agent char. once at the household level 
        IF (ntick.eq.1) THEN
          jinit=1
          do kx=1,8
           do ky=1,20
             n3e=npop2(kx,ky)/4 
		 
c ... here we are setting an agents dthresh which is the amount of time they can last
c ... (10-30 hours) waiting on the road to evac before they give up. 
c ... likewise we set agents departure time (decitiont - i.e. the amt of time they wait before leaving)
c ... failatt is the number of times a vehicle has looked for a spot on road and couldnt find one
c ... destside is either 0 (left side of grid) or 1 (right side of grid)
c ... onerte (for those who enter the southern east-west interstate, indicates which direction they head)
c ... tworte (same thing but for middle east-west interstate connecting tampa and orlando)
c ... exitrow determines which grid cell evacuees exit I75 and I95 

             do n4=1,n3e         
                IF (pevac(kx,ky,n4).ne.1) then 
                   pevac(kx,ky,n4)=0
                end if    
                failatt(kx,ky,n4)=0                   
                dthresh(kx,ky,n4)=9999999
                decisiont(kx,ky,n4)=9999999
                destside(kx,ky,n4)=0
                onerte(kx,ky,n4)=0
                tworte(kx,ky,n4)=0
                exitrow(kx,ky,n4)=0
             end do
           end do 
          end do
c.... counting the spots that are available in upstate Florida (based on tile pop)

          maxspots(1,21)=20000000
          maxspots(8,21)=20000000
          maxspots(2,21)=0
          maxspots(3,21)=0
          maxspots(4,21)=0
          maxspots(5,21)=0
          maxspots(6,21)=0
          maxspots(7,21)=0
          spt(1,21)=0
          spt(4,21)=0
          spt(2,21)=0
          spt(3,21)=0
          spt(5,21)=0
          spt(6,21)=0
          spt(7,21)=0
          spt(8,21)=0
          do lx=1,8
          do ly=8,20
            spt(lx,ly)=0
            dest(lx,ly)=0
            IF (lx.eq.1) then 
               maxspots(lx,ly)=npop(lx,ly)/8
            else if (lx.eq.8) then
               maxspots(lx,ly)=npop(lx,ly)/8
            else if (ly.eq.10) then
               maxspots(lx,ly)=npop(lx,ly)/8  
            else
               maxspots(lx,ly)=0
            end if   
          end do
          end do
        ELSE
          jinit=0
        end if 

c  still more defining initial parameter states		
        if (ntick.eq.2) then
          do kx=1,8
          do ky=1,20
            n3e=npop2(kx,ky)/4      
            do n4=1,n3e
c. mevac=99 means you are eligable for evacuating (i.e., you have a car) 
c. mevac=0			
              if (mevac(kx,ky,n4).eq.99) then
                  if (rand().le.0.01) then
                      mevac(kx,ky,n4)=0
                  else 
                      mevac(kx,ky,n4)=99
                  end if
              end if    
            end do
          end do
          end do 
         call ABM(0,npop2,wrisk,srisk,rrisk,mevac,pevac,
     *     socio,age,nocar,mobl,agewt,
     *     wwt,swt,rwt,barr,mobwt,strmwt,ewt,tstoa)
        end if

        if (mod(ntick,1500).eq.0) then
          do kx=1,8
          do ky=1,20
            n3e=npop2(kx,ky)/4      
c. code that eases agents into evacuation decisions
c. that way there's not an unrealistic big jump in evacuees based on the first forecast

            do n4=1,n3e
              if (mevac(kx,ky,n4).eq.99) then
                if (ntick.le.18000) then 
                  if (rand().le.0.01) then
                      mevac(kx,ky,n4)=0
                  else 
                      mevac(kx,ky,n4)=99
                  end if
                end if
                if (ntick.gt.18000.and.ntick.le.36000) then 
                  if (rand().le.0.02) then
                      mevac(kx,ky,n4)=0
                  else 
                      mevac(kx,ky,n4)=99
                  end if
                end if
                if (ntick.gt.36000.and.ntick.le.54000) then 
                  if (rand().le.0.03) then
                      mevac(kx,ky,n4)=0
                  else 
                      mevac(kx,ky,n4)=99
                  end if
                end if
                if (ntick.gt.54000.and.ntick.le.72000) then 
                  if (rand().le.0.04) then
                      mevac(kx,ky,n4)=0
                  else 
                      mevac(kx,ky,n4)=99
                  end if
                end if
                if (ntick.gt.72000.and.ntick.le.90000) then 
                  if (rand().le.0.05) then
                      mevac(kx,ky,n4)=0
                  else 
                      mevac(kx,ky,n4)=99
                  end if
                end if
                if (ntick.gt.90000.and.ntick.le.108000) then 
                  if (rand().le.0.06) then
                      mevac(kx,ky,n4)=0
                  else 
                      mevac(kx,ky,n4)=99
                  end if
                end if
                if (ntick.gt.108000.and.ntick.le.126000) then 
                  if (rand().le.0.07) then
                      mevac(kx,ky,n4)=0
                  else 
                      mevac(kx,ky,n4)=99
                  end if
                end if
                if (ntick.gt.126000.and.ntick.le.144000) then 
                  if (rand().le.0.08) then
                      mevac(kx,ky,n4)=0
                  else 
                      mevac(kx,ky,n4)=99
                  end if
                end if
                if (ntick.gt.144000) then 
                  if (rand().le.0.09) then
                      mevac(kx,ky,n4)=0
                  else 
                      mevac(kx,ky,n4)=99
                  end if
                end if                
              end if    
            end do
          end do
          end do 
c ... here we call the ABM EVERY 30 MINUTES		  		  
         call ABM(0,npop2,wrisk,srisk,rrisk,mevac,pevac,
     *     socio,age,nocar,mobl,agewt,
     *     wwt,swt,rwt,barr,mobwt,strmwt,ewt,tstoa)
        end if
c   mevac = 1 when the decision is to go. 
c.  mevac = 4 are those that simply cannot leave due to lack of car etc. 		
        do kx=1,8
        do ky=1,20
          n3e=npop2(kx,ky)/4      
          do n4=1,n3e
            if (mevac(kx,ky,n4).eq.1) then
               if (ntick.le.216000) then
                  decisiont(kx,ky,n4)=ntick+int(rand()*42000)
               else if (ntick.gt.216000.and.ntick.le.288000) then
                  decisiont(kx,ky,n4)=ntick+int(rand()*36000)
               else if (ntick.gt.288000.and.ntick.le.360000) then
                  decisiont(kx,ky,n4)=ntick+int(rand()*30000)
               else if (ntick.gt.360000.and.ntick.le.396000) then
                  decisiont(kx,ky,n4)=ntick+int(rand()*24000)
               else
                  decisiont(kx,ky,n4)=ntick+int(rand()*18000)
               end if    
            mevac(kx,ky,n4)=10
            end if
            if (mevac(kx,ky,n4).eq.4.and.evac(kx,ky).eq.1) then
               if (ntick.le.216000) then
                  decisiont(kx,ky,n4)=ntick+int(rand()*42000)
               else if (ntick.gt.216000.and.ntick.le.288000) then
                  decisiont(kx,ky,n4)=ntick+int(rand()*36000)
               else if (ntick.gt.288000.and.ntick.le.360000) then
                  decisiont(kx,ky,n4)=ntick+int(rand()*30000)
               else if (ntick.gt.360000.and.ntick.le.396000) then
                  decisiont(kx,ky,n4)=ntick+int(rand()*24000)
               else
                  decisiont(kx,ky,n4)=ntick+int(rand()*18000)
               end if 
            mevac(kx,ky,n4)=44
            end if            
          end do
        end do 
        end do
         cnt0=0
         cnt1=0
         cnt2=0
         cnt3=0
         cnt4=0
         cnt5=0
         cnt6=0
         cnt7=0
         cnt8=0
         cnt9=0
         cnt99=0
         do kx=1,8
           do ky=1,20
             ct0(kx,ky)=0
             ct1(kx,ky)=0
             ct2(kx,ky)=0
             ct3(kx,ky)=0
             ct4(kx,ky)=0
             ct5(kx,ky)=0
             ct6(kx,ky)=0
             ct7(kx,ky)=0
             ct8(kx,ky)=0  
             ct9(kx,ky)=0
             ct10(kx,ky)=0
             ct11(kx,ky)=0
             ct12(kx,ky)=0
             ct99(kx,ky)=0
             n3e=npop2(kx,ky)/4 
             side=0
             do n4=1,n3e             
c... assigning desintation tiles... and directions... to those wanting to leave 
               if (mevac(kx,ky,n4).eq.10.and. 
     *           ntick.gt.decisiont(kx,ky,n4)) THEN     
                 dthresh(kx,ky,n4)=(1+
     *               int((rand()*tstoa(kx,ky)*1500)))        
                  if (ky.le.2) then
                    if (kx.le.6) then
                         side=1
                    else if (kx.gt.6) then
                      if (rand().le.0.40) then
                         side=1
                      else 
                         side=2
                      end if
                    end if                      
                    if (side.eq.1) then
                      if (rand().le.0.0001) then
                        if (spt(1,21).lt.maxspots(1,21)) then
                           exitrow(kx,ky,n4)=21
                           spt(1,21)=spt(1,21)+1
                           dest(1,21)=spt(1,21)                           
                           onerte(kx,ky,n4)=1
                           tworte(kx,ky,n4)=0
                           goto 623
                        end if
                      else                        

                        do lx=4,1,-1
                        if (spt(lx,10).lt.maxspots(lx,10)) then
                           spt(lx,10)=spt(lx,10)+1
                           dest(lx,10)=spt(lx,10)
                           exitrow(kx,ky,n4)=10
                           onerte(kx,ky,n4)=1
                           tworte(kx,ky,n4)=1
                           goto 623
                        end if
                        end do
		
                        do ly=8,20
                        if (spt(1,ly).lt.maxspots(1,ly)) then
                           spt(1,ly)=spt(1,ly)+1
                           dest(1,ly)=spt(1,ly)
                           exitrow(kx,ky,n4)=ly
                           onerte(kx,ky,n4)=1
                           tworte(kx,ky,n4)=0
                           goto 623
                        end if
                        end do
                      end if  
                    ELSE if (side.eq.2) then
                      if (rand().le.0.0001) then
                        if (spt(8,21).lt.maxspots(8,21)) then
                           exitrow(kx,ky,n4)=21
                           spt(8,21)=spt(8,21)+1
                           dest(8,21)=spt(8,21)                             
                           onerte(kx,ky,n4)=2
                           tworte(kx,ky,n4)=0
                           goto 623
                        end if
                      else 
					  
                        do lx=5,8
                        if (spt(lx,10).lt.maxspots(lx,10)) then
                           spt(lx,10)=spt(lx,10)+1
                           dest(lx,10)=spt(lx,10)
                           exitrow(kx,ky,n4)=10
                           onerte(kx,ky,n4)=2
                           tworte(kx,ky,n4)=2
                           goto 623
                        end if
                        end do
						
                        do ly=8,20
                        if (spt(8,ly).lt.maxspots(8,ly)) then
                           spt(8,ly)=spt(8,ly)+1
                           dest(8,ly)=spt(8,ly)
                           exitrow(kx,ky,n4)=ly
                           onerte(kx,ky,n4)=2
                           tworte(kx,ky,n4)=0
                           goto 623
                        end if
                        end do
                      end if  
                    end if                          
                  ELSE if (ky.ge.3.and.ky.le.7) then
                    IF (kx.le.4) then
                      if (rand().le.0.0001) then
                        if (spt(1,21).lt.maxspots(1,21)) then
                           exitrow(kx,ky,n4)=21
                           spt(1,21)=spt(1,21)+1
                           dest(1,21)=spt(1,21)                             
                           onerte(kx,ky,n4)=0
                           tworte(kx,ky,n4)=0
                           goto 623
                         end if  
                      else   
                        do lx=4,1,-1
                        if (spt(lx,10).lt.maxspots(lx,10)) then
                           spt(lx,10)=spt(lx,10)+1
                           dest(lx,10)=spt(lx,10)
                           exitrow(kx,ky,n4)=10
                           onerte(kx,ky,n4)=0
                           tworte(kx,ky,n4)=1
                           goto 623
                        end if
                        end do
						
                        do ly=8,20
                        if (spt(1,ly).lt.maxspots(1,ly)) then
                           spt(1,ly)=spt(1,ly)+1
                           dest(1,ly)=spt(1,ly)
                           exitrow(kx,ky,n4)=ly
                           onerte(kx,ky,n4)=0
                           tworte(kx,ky,n4)=0
                           goto 623
                        end if
                        end do
                      end if  
                    else
                      if (rand().le.0.00001) then
                        if (spt(8,21).lt.maxspots(8,21)) then
                           exitrow(kx,ky,n4)=21
                           spt(8,21)=spt(8,21)+1
                           dest(8,21)=spt(8,21)                             
                           onerte(kx,ky,n4)=0
                           tworte(kx,ky,n4)=0
                           goto 623
                         end if  
                      else                         
                        do lx=5,8
                        if (spt(lx,10).lt.maxspots(lx,10)) then
                           spt(lx,10)=spt(lx,10)+1
                           dest(lx,10)=spt(lx,10)
                           exitrow(kx,ky,n4)=10
                           onerte(kx,ky,n4)=0
                           tworte(kx,ky,n4)=2
                           goto 623
                        end if
                        end do
                        do ly=8,20
                        if (spt(8,ly).lt.maxspots(8,ly)) then
                           spt(8,ly)=spt(8,ly)+1
                           dest(8,ly)=spt(8,ly)
                           exitrow(kx,ky,n4)=ly
                           onerte(kx,ky,n4)=0
                           tworte(kx,ky,n4)=0
                           goto 623
                        end if
                        end do
                      end if                        
                    end if
                  else 
                    if (kx.le.4) then
                    if (spt(1,21).lt.maxspots(1,21)) then
                           exitrow(kx,ky,n4)=21
                           spt(1,21)=spt(1,21)+1
                           dest(1,21)=spt(1,21)                             
                           onerte(kx,ky,n4)=0
                           tworte(kx,ky,n4)=0
                           goto 623
                    end if 
                    end if 
                    if (kx.gt.4) then
                    if (spt(8,21).lt.maxspots(8,21)) then
                           exitrow(kx,ky,n4)=21
                           spt(8,21)=spt(8,21)+1
                           dest(8,21)=spt(8,21)                             
                           onerte(kx,ky,n4)=0
                           tworte(kx,ky,n4)=0
                           goto 623
                    end if 
                    end if                    
                  end if
               else              
               end if
 623     continue
 
c ... if the decision made in subroutine is GO then check for an open spot on road
c ... if a spot is open then mevac = 2 (agent group is gone) and put a car on the road
               if (ntick.gt.decisiont(kx,ky,n4)) THEN
               if (mevac(kx,ky,n4).eq.10.or.mevac(kx,ky,n4).eq.8) THEN
               ncx=1+int(rand()*real(NCELLX))
               do lx=1,2
               DO le=1,3
               if (ky.eq.1) then			   
               if (onerte(kx,ky,n4).eq.2) then
               if (carse(kx,ky,ncx,le).eq.0.and.
     *           npop(kx,ky).gt.0) THEN
                    mevac(kx,ky,n4)=2
                    carse(kx,ky,ncx,le)=1
                    onese(kx,ky,ncx,le)=onerte(kx,ky,n4)
                    twose(kx,ky,ncx,le)=tworte(kx,ky,n4)
                    exitse(kx,ky,ncx,le)=exitrow(kx,ky,n4)
                    velse(kx,ky,ncx,le)=nVMX
                    identse(kx,ky,ncx,le)=nlast+1
                    npopl=npop(kx,ky)
                    npop(kx,ky)=max(npop(kx,ky)-2,0)
                    num=npopl-npop(kx,ky)
                    carpse(kx,ky,ncx,le)=num
                    nlast=nlast+1    
               end if
               end if
               if (onerte(kx,ky,n4).eq.1) then
               if (carsw(kx,ky,ncx,le).eq.0.and.
     *           npop(kx,ky).gt.0) THEN
                    mevac(kx,ky,n4)=2
                    carsw(kx,ky,ncx,le)=1
                    onesw(kx,ky,ncx,le)=onerte(kx,ky,n4)
                    twosw(kx,ky,ncx,le)=tworte(kx,ky,n4)
                    exitsw(kx,ky,ncx,le)=exitrow(kx,ky,n4)
                    velsw(kx,ky,ncx,le)=nVMX
                    identsw(kx,ky,ncx,le)=nlast+1
                    npopl=npop(kx,ky)
                    npop(kx,ky)=max(npop(kx,ky)-2,0)
                    num=npopl-npop(kx,ky)
                    carpsw(kx,ky,ncx,le)=num
                    nlast=nlast+1                
               end if
               end if                
			   end if 
			   
               if (ky.eq.2) then			   
               if (onerte(kx,ky,n4).eq.2) then
               if (carse(kx,ky-1,ncx,le).eq.0.and.
     *           npop(kx,ky).gt.0) THEN
                    mevac(kx,ky,n4)=2
                    carse(kx,ky-1,ncx,le)=1
                    onese(kx,ky-1,ncx,le)=onerte(kx,ky,n4)
                    twose(kx,ky-1,ncx,le)=tworte(kx,ky,n4)
                    exitse(kx,ky-1,ncx,le)=exitrow(kx,ky,n4)
                    velse(kx,ky-1,ncx,le)=nVMX
                    identse(kx,ky-1,ncx,le)=nlast+1
                    npopl=npop(kx,ky)
                    npop(kx,ky)=max(npop(kx,ky)-2,0)
                    num=npopl-npop(kx,ky)
                    carpse(kx,ky-1,ncx,le)=num
                    nlast=nlast+1    
               end if
               end if
               if (onerte(kx,ky,n4).eq.1) then
               if (carsw(kx,ky-1,ncx,le).eq.0.and.
     *           npop(kx,ky).gt.0) THEN
                    mevac(kx,ky,n4)=2
                    carsw(kx,ky-1,ncx,le)=1
                    onesw(kx,ky-1,ncx,le)=onerte(kx,ky,n4)
                    twosw(kx,ky-1,ncx,le)=tworte(kx,ky,n4)
                    exitsw(kx,ky-1,ncx,le)=exitrow(kx,ky,n4)
                    velsw(kx,ky-1,ncx,le)=nVMX
                    identsw(kx,ky-1,ncx,le)=nlast+1
                    npopl=npop(kx,ky)
                    npop(kx,ky)=max(npop(kx,ky)-2,0)
                    num=npopl-npop(kx,ky)
                    carpsw(kx,ky,ncx,le)=num
                    nlast=nlast+1                
               end if
               end if                
			   end if			   
			   
               if (ky.eq.4) then
				if (onerte(kx,ky,n4).eq.0) then
				if (care(kx,ky-1,ncx,lx).eq.0.and.
     *           	npop(kx,ky).gt.0) THEN
					mevac(kx,ky,n4)=2
					care(kx,ky-1,ncx,lx)=1
					oneh(kx,ky-1,ncx,lx)=onerte(kx,ky,n4)
					twoh(kx,ky-1,ncx,lx)=tworte(kx,ky,n4)
					exith(kx,ky-1,ncx,lx)=exitrow(kx,ky,n4)
					vele(kx,ky-1,ncx,lx)=nVMX2
					idente(kx,ky-1,ncx,lx)=nlast+1
					npopl=npop(kx,ky)
					npop(kx,ky)=max(npop(kx,ky)-2,0)
					num=npopl-npop(kx,ky)
					carpe(kx,ky-1,ncx,lx)=num
					nlast=nlast+1
				end if 
				end if	   
			   else if (ky.eq.6) then		   
				if (onerte(kx,ky,n4).eq.0) then
				if (care(kx,ky-1,ncx,lx).eq.0.and.
     *           	npop(kx,ky).gt.0) THEN
					mevac(kx,ky,n4)=2
					care(kx,ky-1,ncx,lx)=1
					oneh(kx,ky-1,ncx,lx)=onerte(kx,ky,n4)
					twoh(kx,ky-1,ncx,lx)=tworte(kx,ky,n4)
					exith(kx,ky-1,ncx,lx)=exitrow(kx,ky,n4)
					vele(kx,ky-1,ncx,lx)=nVMX2
					idente(kx,ky-1,ncx,lx)=nlast+1
					npopl=npop(kx,ky)
					npop(kx,ky)=max(npop(kx,ky)-2,0)
					num=npopl-npop(kx,ky)
					carpe(kx,ky-1,ncx,lx)=num
					nlast=nlast+1
				end if 
				end if				   
			   else if (ky.eq.8) then		   
				if (onerte(kx,ky,n4).eq.0) then
				if (care(kx,ky-1,ncx,lx).eq.0.and.
     *           	npop(kx,ky).gt.0) THEN
					mevac(kx,ky,n4)=2
					care(kx,ky-1,ncx,lx)=1
					oneh(kx,ky-1,ncx,lx)=onerte(kx,ky,n4)
					twoh(kx,ky-1,ncx,lx)=tworte(kx,ky,n4)
					exith(kx,ky-1,ncx,lx)=exitrow(kx,ky,n4)
					vele(kx,ky-1,ncx,lx)=nVMX2
					idente(kx,ky-1,ncx,lx)=nlast+1
					npopl=npop(kx,ky)
					npop(kx,ky)=max(npop(kx,ky)-2,0)
					num=npopl-npop(kx,ky)
					carpe(kx,ky-1,ncx,lx)=num
					nlast=nlast+1
				end if 
				end if	
			   else if (ky.eq.10) then		   
				if (onerte(kx,ky,n4).eq.0) then
				if (care(kx,ky-1,ncx,lx).eq.0.and.
     *           	npop(kx,ky).gt.0) THEN
					mevac(kx,ky,n4)=2
					care(kx,ky-1,ncx,lx)=1
					oneh(kx,ky-1,ncx,lx)=onerte(kx,ky,n4)
					twoh(kx,ky-1,ncx,lx)=tworte(kx,ky,n4)
					exith(kx,ky-1,ncx,lx)=exitrow(kx,ky,n4)
					vele(kx,ky-1,ncx,lx)=nVMX2
					idente(kx,ky-1,ncx,lx)=nlast+1
					npopl=npop(kx,ky)
					npop(kx,ky)=max(npop(kx,ky)-2,0)
					num=npopl-npop(kx,ky)
					carpe(kx,ky-1,ncx,lx)=num
					nlast=nlast+1
				end if 
				end if	
			   else if (ky.eq.12) then		   
				if (onerte(kx,ky,n4).eq.0) then
				if (care(kx,ky-1,ncx,lx).eq.0.and.
     *           	npop(kx,ky).gt.0) THEN
					mevac(kx,ky,n4)=2
					care(kx,ky-1,ncx,lx)=1
					oneh(kx,ky-1,ncx,lx)=onerte(kx,ky,n4)
					twoh(kx,ky-1,ncx,lx)=tworte(kx,ky,n4)
					exith(kx,ky-1,ncx,lx)=exitrow(kx,ky,n4)
					vele(kx,ky-1,ncx,lx)=nVMX2
					idente(kx,ky-1,ncx,lx)=nlast+1
					npopl=npop(kx,ky)
					npop(kx,ky)=max(npop(kx,ky)-2,0)
					num=npopl-npop(kx,ky)
					carpe(kx,ky-1,ncx,lx)=num
					nlast=nlast+1
				end if 
				end if	
			   else if (ky.eq.14) then		   
				if (onerte(kx,ky,n4).eq.0) then
				if (care(kx,ky-1,ncx,lx).eq.0.and.
     *           	npop(kx,ky).gt.0) THEN
					mevac(kx,ky,n4)=2
					care(kx,ky-1,ncx,lx)=1
					oneh(kx,ky-1,ncx,lx)=onerte(kx,ky,n4)
					twoh(kx,ky-1,ncx,lx)=tworte(kx,ky,n4)
					exith(kx,ky-1,ncx,lx)=exitrow(kx,ky,n4)
					vele(kx,ky-1,ncx,lx)=nVMX2
					idente(kx,ky-1,ncx,lx)=nlast+1
					npopl=npop(kx,ky)
					npop(kx,ky)=max(npop(kx,ky)-2,0)
					num=npopl-npop(kx,ky)
					carpe(kx,ky-1,ncx,lx)=num
					nlast=nlast+1
				end if 
				end if	
			   else if (ky.eq.16) then		   
				if (onerte(kx,ky,n4).eq.0) then
				if (care(kx,ky-1,ncx,lx).eq.0.and.
     *           	npop(kx,ky).gt.0) THEN
					mevac(kx,ky,n4)=2
					care(kx,ky-1,ncx,lx)=1
					oneh(kx,ky-1,ncx,lx)=onerte(kx,ky,n4)
					twoh(kx,ky-1,ncx,lx)=tworte(kx,ky,n4)
					exith(kx,ky-1,ncx,lx)=exitrow(kx,ky,n4)
					vele(kx,ky-1,ncx,lx)=nVMX2
					idente(kx,ky-1,ncx,lx)=nlast+1
					npopl=npop(kx,ky)
					npop(kx,ky)=max(npop(kx,ky)-2,0)
					num=npopl-npop(kx,ky)
					carpe(kx,ky-1,ncx,lx)=num
					nlast=nlast+1
				end if 
				end if		
			   else if (ky.eq.18) then		   
				if (onerte(kx,ky,n4).eq.0) then
				if (care(kx,ky-1,ncx,lx).eq.0.and.
     *           	npop(kx,ky).gt.0) THEN
					mevac(kx,ky,n4)=2
					care(kx,ky-1,ncx,lx)=1
					oneh(kx,ky-1,ncx,lx)=onerte(kx,ky,n4)
					twoh(kx,ky-1,ncx,lx)=tworte(kx,ky,n4)
					exith(kx,ky-1,ncx,lx)=exitrow(kx,ky,n4)
					vele(kx,ky-1,ncx,lx)=nVMX2
					idente(kx,ky-1,ncx,lx)=nlast+1
					npopl=npop(kx,ky)
					npop(kx,ky)=max(npop(kx,ky)-2,0)
					num=npopl-npop(kx,ky)
					carpe(kx,ky-1,ncx,lx)=num
					nlast=nlast+1
				end if 
				end if
			   else if (ky.eq.20) then		   
				if (onerte(kx,ky,n4).eq.0) then
				if (care(kx,ky-1,ncx,lx).eq.0.and.
     *           	npop(kx,ky).gt.0) THEN
					mevac(kx,ky,n4)=2
					care(kx,ky-1,ncx,lx)=1
					oneh(kx,ky-1,ncx,lx)=onerte(kx,ky,n4)
					twoh(kx,ky-1,ncx,lx)=tworte(kx,ky,n4)
					exith(kx,ky-1,ncx,lx)=exitrow(kx,ky,n4)
					vele(kx,ky-1,ncx,lx)=nVMX2
					idente(kx,ky-1,ncx,lx)=nlast+1
					npopl=npop(kx,ky)
					npop(kx,ky)=max(npop(kx,ky)-2,0)
					num=npopl-npop(kx,ky)
					carpe(kx,ky-1,ncx,lx)=num
					nlast=nlast+1
				end if 
				end if				
			   else 
				if (onerte(kx,ky,n4).eq.0) then
				if (care(kx,ky,ncx,lx).eq.0.and.
     *           	npop(kx,ky).gt.0) THEN
					mevac(kx,ky,n4)=2
					care(kx,ky,ncx,lx)=1
					oneh(kx,ky,ncx,lx)=onerte(kx,ky,n4)
					twoh(kx,ky,ncx,lx)=tworte(kx,ky,n4)
					exith(kx,ky,ncx,lx)=exitrow(kx,ky,n4)
					vele(kx,ky,ncx,lx)=nVMX2
					idente(kx,ky,ncx,lx)=nlast+1
					npopl=npop(kx,ky)
					npop(kx,ky)=max(npop(kx,ky)-2,0)
					num=npopl-npop(kx,ky)
					carpe(kx,ky,ncx,lx)=num
					nlast=nlast+1
				end if 
				end if	   
               end if			   
               end do
			   end do
               failatt(kx,ky,n4)=failatt(kx,ky,n4)+1
               end if
               end if
c..... if mevac=1 at dec time to become mevac(kx,ky,n4)=2, then mevac(kx,ky,n4)=8
               if (mevac(kx,ky,n4).eq.10.and. 
     *           ntick.gt.decisiont(kx,ky,n4)) THEN
                  mevac(kx,ky,n4)=8
               END IF   
c ... if that car couldn't find a spot open, then count the number of ticks they're waiting               
c ... if they're wait time matches their max patience, then they'll shelter at home
               if (failatt(kx,ky,n4).eq.dthresh(kx,ky,n4)) then 
                 mevac(kx,ky,n4)=3
               end if

               if (mevac(kx,ky,n4).eq.44.and. 
     *           ntick.gt.decisiont(kx,ky,n4)) THEN
                  mevac(kx,ky,n4)=5
               END IF   

c ... keeping a talley on agent decisions, will print to screen below               
               if (mevac(kx,ky,n4).eq.0) then 
                  cnt0=cnt0+1
                  IF (ky.eq.1.and.kx.eq.1) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.1.and.kx.eq.2) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.1.and.kx.eq.3) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.1.and.kx.eq.4) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.1.and.kx.eq.5) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.1.and.kx.eq.6) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.1.and.kx.eq.7) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.1.and.kx.eq.8) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.2.and.kx.eq.1) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.2.and.kx.eq.2) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.2.and.kx.eq.3) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.2.and.kx.eq.4) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.2.and.kx.eq.5) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.2.and.kx.eq.6) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.2.and.kx.eq.7) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.2.and.kx.eq.8) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.3.and.kx.eq.1) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.3.and.kx.eq.2) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.3.and.kx.eq.3) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.3.and.kx.eq.4) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.3.and.kx.eq.5) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.3.and.kx.eq.6) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.3.and.kx.eq.7) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.3.and.kx.eq.8) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.4.and.kx.eq.1) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.4.and.kx.eq.2) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.4.and.kx.eq.3) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.4.and.kx.eq.4) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.4.and.kx.eq.5) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.4.and.kx.eq.6) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.4.and.kx.eq.7) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.4.and.kx.eq.8) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.5.and.kx.eq.1) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.5.and.kx.eq.2) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.5.and.kx.eq.3) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.5.and.kx.eq.4) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.5.and.kx.eq.5) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.5.and.kx.eq.6) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.5.and.kx.eq.7) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.5.and.kx.eq.8) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.6.and.kx.eq.1) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.6.and.kx.eq.2) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.6.and.kx.eq.3) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.6.and.kx.eq.4) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.6.and.kx.eq.5) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.6.and.kx.eq.6) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.6.and.kx.eq.7) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.6.and.kx.eq.8) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.7.and.kx.eq.1) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.7.and.kx.eq.2) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.7.and.kx.eq.3) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.7.and.kx.eq.4) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.7.and.kx.eq.5) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.7.and.kx.eq.6) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.7.and.kx.eq.7) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.7.and.kx.eq.8) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.8.and.kx.eq.1) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.8.and.kx.eq.2) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.8.and.kx.eq.3) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.8.and.kx.eq.4) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.8.and.kx.eq.5) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.8.and.kx.eq.6) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.8.and.kx.eq.7) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.8.and.kx.eq.8) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.9.and.kx.eq.1) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.9.and.kx.eq.2) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.9.and.kx.eq.3) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.9.and.kx.eq.4) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.9.and.kx.eq.5) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.9.and.kx.eq.6) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.9.and.kx.eq.7) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.9.and.kx.eq.8) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.10.and.kx.eq.1) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.10.and.kx.eq.2) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.10.and.kx.eq.3) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.10.and.kx.eq.4) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.10.and.kx.eq.5) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.10.and.kx.eq.6) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.10.and.kx.eq.7) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.10.and.kx.eq.8) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else if (ky.eq.11.and.kx.eq.1) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.11.and.kx.eq.2) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.11.and.kx.eq.3) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.11.and.kx.eq.4) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.11.and.kx.eq.5) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.11.and.kx.eq.6) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.11.and.kx.eq.7) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.11.and.kx.eq.8) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.12.and.kx.eq.1) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.12.and.kx.eq.2) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.12.and.kx.eq.3) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.12.and.kx.eq.4) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.12.and.kx.eq.5) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.12.and.kx.eq.6) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.12.and.kx.eq.7) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.12.and.kx.eq.8) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.13.and.kx.eq.1) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.13.and.kx.eq.2) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.13.and.kx.eq.3) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.13.and.kx.eq.4) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.13.and.kx.eq.5) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.13.and.kx.eq.6) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.13.and.kx.eq.7) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.13.and.kx.eq.8) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.14.and.kx.eq.1) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.14.and.kx.eq.2) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.14.and.kx.eq.3) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.14.and.kx.eq.4) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.14.and.kx.eq.5) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.14.and.kx.eq.6) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.14.and.kx.eq.7) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.14.and.kx.eq.8) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.15.and.kx.eq.1) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.15.and.kx.eq.2) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.15.and.kx.eq.3) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.15.and.kx.eq.4) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.15.and.kx.eq.5) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.15.and.kx.eq.6) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.15.and.kx.eq.7) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.15.and.kx.eq.8) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.16.and.kx.eq.1) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.16.and.kx.eq.2) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.16.and.kx.eq.3) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.16.and.kx.eq.4) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.16.and.kx.eq.5) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.16.and.kx.eq.6) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.16.and.kx.eq.7) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.16.and.kx.eq.8) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.17.and.kx.eq.1) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.17.and.kx.eq.2) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.17.and.kx.eq.3) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.17.and.kx.eq.4) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.17.and.kx.eq.5) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.17.and.kx.eq.6) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.17.and.kx.eq.7) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.17.and.kx.eq.8) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.18.and.kx.eq.1) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.18.and.kx.eq.2) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.18.and.kx.eq.3) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.18.and.kx.eq.4) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.18.and.kx.eq.5) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.18.and.kx.eq.6) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.18.and.kx.eq.7) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.18.and.kx.eq.8) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.19.and.kx.eq.1) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.19.and.kx.eq.2) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.19.and.kx.eq.3) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.19.and.kx.eq.4) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.19.and.kx.eq.5) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.19.and.kx.eq.6) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.19.and.kx.eq.7) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.19.and.kx.eq.8) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.20.and.kx.eq.1) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.20.and.kx.eq.2) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.20.and.kx.eq.3) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.20.and.kx.eq.4) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.20.and.kx.eq.5) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.20.and.kx.eq.6) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.20.and.kx.eq.7) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  else IF (ky.eq.20.and.kx.eq.8) then
                    ct0(kx,ky)=ct0(kx,ky)+1
                  end if 
               end if 
               if (mevac(kx,ky,n4).eq.1.or.mevac(kx,ky,n4).eq.10) then 
                  cnt1=cnt1+1
                  IF (ky.eq.1.and.kx.eq.1) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.1.and.kx.eq.2) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.1.and.kx.eq.3) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.1.and.kx.eq.4) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.1.and.kx.eq.5) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.1.and.kx.eq.6) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.1.and.kx.eq.7) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.1.and.kx.eq.8) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.2.and.kx.eq.1) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.2.and.kx.eq.2) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.2.and.kx.eq.3) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.2.and.kx.eq.4) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.2.and.kx.eq.5) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.2.and.kx.eq.6) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.2.and.kx.eq.7) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.2.and.kx.eq.8) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.3.and.kx.eq.1) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.3.and.kx.eq.2) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.3.and.kx.eq.3) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.3.and.kx.eq.4) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.3.and.kx.eq.5) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.3.and.kx.eq.6) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.3.and.kx.eq.7) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.3.and.kx.eq.8) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.4.and.kx.eq.1) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.4.and.kx.eq.2) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.4.and.kx.eq.3) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.4.and.kx.eq.4) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.4.and.kx.eq.5) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.4.and.kx.eq.6) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.4.and.kx.eq.7) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.4.and.kx.eq.8) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.5.and.kx.eq.1) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.5.and.kx.eq.2) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.5.and.kx.eq.3) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.5.and.kx.eq.4) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.5.and.kx.eq.5) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.5.and.kx.eq.6) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.5.and.kx.eq.7) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.5.and.kx.eq.8) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.6.and.kx.eq.1) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.6.and.kx.eq.2) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.6.and.kx.eq.3) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.6.and.kx.eq.4) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.6.and.kx.eq.5) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.6.and.kx.eq.6) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.6.and.kx.eq.7) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.6.and.kx.eq.8) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.7.and.kx.eq.1) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.7.and.kx.eq.2) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.7.and.kx.eq.3) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.7.and.kx.eq.4) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.7.and.kx.eq.5) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.7.and.kx.eq.6) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.7.and.kx.eq.7) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.7.and.kx.eq.8) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.8.and.kx.eq.1) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.8.and.kx.eq.2) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.8.and.kx.eq.3) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.8.and.kx.eq.4) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.8.and.kx.eq.5) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.8.and.kx.eq.6) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.8.and.kx.eq.7) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.8.and.kx.eq.8) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.9.and.kx.eq.1) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.9.and.kx.eq.2) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.9.and.kx.eq.3) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.9.and.kx.eq.4) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.9.and.kx.eq.5) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.9.and.kx.eq.6) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.9.and.kx.eq.7) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.9.and.kx.eq.8) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.10.and.kx.eq.1) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.10.and.kx.eq.2) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.10.and.kx.eq.3) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.10.and.kx.eq.4) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.10.and.kx.eq.5) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.10.and.kx.eq.6) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.10.and.kx.eq.7) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.10.and.kx.eq.8) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.11.and.kx.eq.1) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.11.and.kx.eq.2) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.11.and.kx.eq.3) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.11.and.kx.eq.4) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.11.and.kx.eq.5) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.11.and.kx.eq.6) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.11.and.kx.eq.7) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.11.and.kx.eq.8) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.12.and.kx.eq.1) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.12.and.kx.eq.2) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.12.and.kx.eq.3) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.12.and.kx.eq.4) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.12.and.kx.eq.5) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.12.and.kx.eq.6) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.12.and.kx.eq.7) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.12.and.kx.eq.8) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.13.and.kx.eq.1) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.13.and.kx.eq.2) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.13.and.kx.eq.3) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.13.and.kx.eq.4) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.13.and.kx.eq.5) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.13.and.kx.eq.6) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.13.and.kx.eq.7) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.13.and.kx.eq.8) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.14.and.kx.eq.1) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.14.and.kx.eq.2) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.14.and.kx.eq.3) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.14.and.kx.eq.4) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.14.and.kx.eq.5) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.14.and.kx.eq.6) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.14.and.kx.eq.7) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.14.and.kx.eq.8) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.15.and.kx.eq.1) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.15.and.kx.eq.2) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.15.and.kx.eq.3) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.15.and.kx.eq.4) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.15.and.kx.eq.5) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.15.and.kx.eq.6) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.15.and.kx.eq.7) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.15.and.kx.eq.8) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.16.and.kx.eq.1) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.16.and.kx.eq.2) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.16.and.kx.eq.3) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.16.and.kx.eq.4) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.16.and.kx.eq.5) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.16.and.kx.eq.6) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.16.and.kx.eq.7) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.16.and.kx.eq.8) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.17.and.kx.eq.1) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.17.and.kx.eq.2) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.17.and.kx.eq.3) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.17.and.kx.eq.4) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.17.and.kx.eq.5) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.17.and.kx.eq.6) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.17.and.kx.eq.7) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.17.and.kx.eq.8) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.18.and.kx.eq.1) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.18.and.kx.eq.2) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.18.and.kx.eq.3) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.18.and.kx.eq.4) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.18.and.kx.eq.5) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.18.and.kx.eq.6) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.18.and.kx.eq.7) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.18.and.kx.eq.8) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.19.and.kx.eq.1) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.19.and.kx.eq.2) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.19.and.kx.eq.3) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.19.and.kx.eq.4) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.19.and.kx.eq.5) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.19.and.kx.eq.6) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.19.and.kx.eq.7) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.19.and.kx.eq.8) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.20.and.kx.eq.1) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.20.and.kx.eq.2) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.20.and.kx.eq.3) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.20.and.kx.eq.4) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.20.and.kx.eq.5) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.20.and.kx.eq.6) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.20.and.kx.eq.7) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  else IF (ky.eq.20.and.kx.eq.8) then
                    ct1(kx,ky)=ct1(kx,ky)+1
                  end if
               end if 
               if (mevac(kx,ky,n4).eq.2) then 
                  cnt2=cnt2+1  
                  IF (ky.eq.1.and.kx.eq.1) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.1.and.kx.eq.2) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.1.and.kx.eq.3) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.1.and.kx.eq.4) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.1.and.kx.eq.5) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.1.and.kx.eq.6) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.1.and.kx.eq.7) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.1.and.kx.eq.8) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.2.and.kx.eq.1) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.2.and.kx.eq.2) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.2.and.kx.eq.3) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.2.and.kx.eq.4) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.2.and.kx.eq.5) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.2.and.kx.eq.6) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.2.and.kx.eq.7) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.2.and.kx.eq.8) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.3.and.kx.eq.1) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.3.and.kx.eq.2) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.3.and.kx.eq.3) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.3.and.kx.eq.4) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.3.and.kx.eq.5) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.3.and.kx.eq.6) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.3.and.kx.eq.7) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.3.and.kx.eq.8) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.4.and.kx.eq.1) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.4.and.kx.eq.2) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.4.and.kx.eq.3) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.4.and.kx.eq.4) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.4.and.kx.eq.5) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.4.and.kx.eq.6) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.4.and.kx.eq.7) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.4.and.kx.eq.8) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.5.and.kx.eq.1) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.5.and.kx.eq.2) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.5.and.kx.eq.3) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.5.and.kx.eq.4) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.5.and.kx.eq.5) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.5.and.kx.eq.6) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.5.and.kx.eq.7) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.5.and.kx.eq.8) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.6.and.kx.eq.1) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.6.and.kx.eq.2) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.6.and.kx.eq.3) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.6.and.kx.eq.4) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.6.and.kx.eq.5) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.6.and.kx.eq.6) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.6.and.kx.eq.7) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.6.and.kx.eq.8) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.7.and.kx.eq.1) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.7.and.kx.eq.2) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.7.and.kx.eq.3) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.7.and.kx.eq.4) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.7.and.kx.eq.5) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.7.and.kx.eq.6) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.7.and.kx.eq.7) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.7.and.kx.eq.8) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.8.and.kx.eq.1) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.8.and.kx.eq.2) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.8.and.kx.eq.3) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.8.and.kx.eq.4) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.8.and.kx.eq.5) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.8.and.kx.eq.6) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.8.and.kx.eq.7) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.8.and.kx.eq.8) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.9.and.kx.eq.1) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.9.and.kx.eq.2) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.9.and.kx.eq.3) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.9.and.kx.eq.4) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.9.and.kx.eq.5) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.9.and.kx.eq.6) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.9.and.kx.eq.7) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.9.and.kx.eq.8) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.10.and.kx.eq.1) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.10.and.kx.eq.2) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.10.and.kx.eq.3) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.10.and.kx.eq.4) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.10.and.kx.eq.5) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.10.and.kx.eq.6) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.10.and.kx.eq.7) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.10.and.kx.eq.8) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.11.and.kx.eq.1) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.11.and.kx.eq.2) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.11.and.kx.eq.3) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.11.and.kx.eq.4) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.11.and.kx.eq.5) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.11.and.kx.eq.6) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.11.and.kx.eq.7) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.11.and.kx.eq.8) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.12.and.kx.eq.1) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.12.and.kx.eq.2) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.12.and.kx.eq.3) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.12.and.kx.eq.4) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.12.and.kx.eq.5) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.12.and.kx.eq.6) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.12.and.kx.eq.7) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.12.and.kx.eq.8) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.13.and.kx.eq.1) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.13.and.kx.eq.2) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.13.and.kx.eq.3) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.13.and.kx.eq.4) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.13.and.kx.eq.5) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.13.and.kx.eq.6) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.13.and.kx.eq.7) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.13.and.kx.eq.8) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.14.and.kx.eq.1) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.14.and.kx.eq.2) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.14.and.kx.eq.3) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.14.and.kx.eq.4) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.14.and.kx.eq.5) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.14.and.kx.eq.6) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.14.and.kx.eq.7) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.14.and.kx.eq.8) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.15.and.kx.eq.1) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.15.and.kx.eq.2) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.15.and.kx.eq.3) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.15.and.kx.eq.4) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.15.and.kx.eq.5) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.15.and.kx.eq.6) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.15.and.kx.eq.7) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.15.and.kx.eq.8) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.16.and.kx.eq.1) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.16.and.kx.eq.2) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.16.and.kx.eq.3) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.16.and.kx.eq.4) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.16.and.kx.eq.5) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.16.and.kx.eq.6) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.16.and.kx.eq.7) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.16.and.kx.eq.8) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.17.and.kx.eq.1) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.17.and.kx.eq.2) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.17.and.kx.eq.3) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.17.and.kx.eq.4) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.17.and.kx.eq.5) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.17.and.kx.eq.6) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.17.and.kx.eq.7) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.17.and.kx.eq.8) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.18.and.kx.eq.1) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.18.and.kx.eq.2) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.18.and.kx.eq.3) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.18.and.kx.eq.4) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.18.and.kx.eq.5) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.18.and.kx.eq.6) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.18.and.kx.eq.7) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.18.and.kx.eq.8) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.19.and.kx.eq.1) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.19.and.kx.eq.2) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.19.and.kx.eq.3) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.19.and.kx.eq.4) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.19.and.kx.eq.5) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.19.and.kx.eq.6) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.19.and.kx.eq.7) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.19.and.kx.eq.8) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.20.and.kx.eq.1) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.20.and.kx.eq.2) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.20.and.kx.eq.3) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.20.and.kx.eq.4) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.20.and.kx.eq.5) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.20.and.kx.eq.6) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.20.and.kx.eq.7) then
                    ct2(kx,ky)=ct2(kx,ky)+1
                  else IF (ky.eq.20.and.kx.eq.8) then
                    ct2(kx,ky)=ct2(kx,ky)+1 
		end if  
               end if 
               if (mevac(kx,ky,n4).eq.3) then 
                  cnt3=cnt3+1  
                  IF (ky.eq.1.and.kx.eq.1) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.1.and.kx.eq.2) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.1.and.kx.eq.3) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.1.and.kx.eq.4) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.1.and.kx.eq.5) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.1.and.kx.eq.6) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.1.and.kx.eq.7) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.1.and.kx.eq.8) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.2.and.kx.eq.1) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.2.and.kx.eq.2) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.2.and.kx.eq.3) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.2.and.kx.eq.4) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.2.and.kx.eq.5) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.2.and.kx.eq.6) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.2.and.kx.eq.7) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.2.and.kx.eq.8) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.3.and.kx.eq.1) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.3.and.kx.eq.2) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.3.and.kx.eq.3) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.3.and.kx.eq.4) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.3.and.kx.eq.5) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.3.and.kx.eq.6) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.3.and.kx.eq.7) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.3.and.kx.eq.8) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.4.and.kx.eq.1) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.4.and.kx.eq.2) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.4.and.kx.eq.3) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.4.and.kx.eq.4) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.4.and.kx.eq.5) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.4.and.kx.eq.6) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.4.and.kx.eq.7) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.4.and.kx.eq.8) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.5.and.kx.eq.1) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.5.and.kx.eq.2) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.5.and.kx.eq.3) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.5.and.kx.eq.4) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.5.and.kx.eq.5) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.5.and.kx.eq.6) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.5.and.kx.eq.7) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.5.and.kx.eq.8) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.6.and.kx.eq.1) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.6.and.kx.eq.2) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.6.and.kx.eq.3) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.6.and.kx.eq.4) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.6.and.kx.eq.5) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.6.and.kx.eq.6) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.6.and.kx.eq.7) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.6.and.kx.eq.8) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.7.and.kx.eq.1) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.7.and.kx.eq.2) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.7.and.kx.eq.3) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.7.and.kx.eq.4) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.7.and.kx.eq.5) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.7.and.kx.eq.6) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.7.and.kx.eq.7) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.7.and.kx.eq.8) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.8.and.kx.eq.1) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.8.and.kx.eq.2) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.8.and.kx.eq.3) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.8.and.kx.eq.4) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.8.and.kx.eq.5) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.8.and.kx.eq.6) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.8.and.kx.eq.7) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.8.and.kx.eq.8) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.9.and.kx.eq.1) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.9.and.kx.eq.2) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.9.and.kx.eq.3) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.9.and.kx.eq.4) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.9.and.kx.eq.5) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.9.and.kx.eq.6) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.9.and.kx.eq.7) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.9.and.kx.eq.8) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.10.and.kx.eq.1) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.10.and.kx.eq.2) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.10.and.kx.eq.3) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.10.and.kx.eq.4) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.10.and.kx.eq.5) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.10.and.kx.eq.6) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.10.and.kx.eq.7) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.10.and.kx.eq.8) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.11.and.kx.eq.1) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.11.and.kx.eq.2) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.11.and.kx.eq.3) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.11.and.kx.eq.4) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.11.and.kx.eq.5) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.11.and.kx.eq.6) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.11.and.kx.eq.7) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.11.and.kx.eq.8) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.12.and.kx.eq.1) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.12.and.kx.eq.2) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.12.and.kx.eq.3) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.12.and.kx.eq.4) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.12.and.kx.eq.5) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.12.and.kx.eq.6) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.12.and.kx.eq.7) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.12.and.kx.eq.8) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.13.and.kx.eq.1) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.13.and.kx.eq.2) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.13.and.kx.eq.3) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.13.and.kx.eq.4) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.13.and.kx.eq.5) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.13.and.kx.eq.6) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.13.and.kx.eq.7) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.13.and.kx.eq.8) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.14.and.kx.eq.1) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.14.and.kx.eq.2) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.14.and.kx.eq.3) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.14.and.kx.eq.4) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.14.and.kx.eq.5) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.14.and.kx.eq.6) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.14.and.kx.eq.7) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.14.and.kx.eq.8) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.15.and.kx.eq.1) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.15.and.kx.eq.2) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.15.and.kx.eq.3) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.15.and.kx.eq.4) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.15.and.kx.eq.5) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.15.and.kx.eq.6) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.15.and.kx.eq.7) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.15.and.kx.eq.8) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.16.and.kx.eq.1) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.16.and.kx.eq.2) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.16.and.kx.eq.3) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.16.and.kx.eq.4) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.16.and.kx.eq.5) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.16.and.kx.eq.6) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.16.and.kx.eq.7) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.16.and.kx.eq.8) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.17.and.kx.eq.1) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.17.and.kx.eq.2) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.17.and.kx.eq.3) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.17.and.kx.eq.4) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.17.and.kx.eq.5) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.17.and.kx.eq.6) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.17.and.kx.eq.7) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.17.and.kx.eq.8) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.18.and.kx.eq.1) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.18.and.kx.eq.2) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.18.and.kx.eq.3) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.18.and.kx.eq.4) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.18.and.kx.eq.5) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.18.and.kx.eq.6) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.18.and.kx.eq.7) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.18.and.kx.eq.8) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.19.and.kx.eq.1) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.19.and.kx.eq.2) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.19.and.kx.eq.3) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.19.and.kx.eq.4) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.19.and.kx.eq.5) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.19.and.kx.eq.6) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.19.and.kx.eq.7) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.19.and.kx.eq.8) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.20.and.kx.eq.1) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.20.and.kx.eq.2) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.20.and.kx.eq.3) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.20.and.kx.eq.4) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.20.and.kx.eq.5) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.20.and.kx.eq.6) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.20.and.kx.eq.7) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                  else IF (ky.eq.20.and.kx.eq.8) then
                    ct3(kx,ky)=ct3(kx,ky)+1
                 end if                            
               end if
               if (mevac(kx,ky,n4).eq.4.or.mevac(kx,ky,n4).eq.44) then 
                  cnt4=cnt4+1 
                  IF (ky.eq.1.and.kx.eq.1) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.1.and.kx.eq.2) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.1.and.kx.eq.3) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.1.and.kx.eq.4) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.1.and.kx.eq.5) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.1.and.kx.eq.6) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.1.and.kx.eq.7) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.1.and.kx.eq.8) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.2.and.kx.eq.1) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.2.and.kx.eq.2) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.2.and.kx.eq.3) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.2.and.kx.eq.4) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.2.and.kx.eq.5) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.2.and.kx.eq.6) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.2.and.kx.eq.7) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.2.and.kx.eq.8) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.3.and.kx.eq.1) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.3.and.kx.eq.2) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.3.and.kx.eq.3) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.3.and.kx.eq.4) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.3.and.kx.eq.5) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.3.and.kx.eq.6) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.3.and.kx.eq.7) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.3.and.kx.eq.8) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.4.and.kx.eq.1) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.4.and.kx.eq.2) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.4.and.kx.eq.3) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.4.and.kx.eq.4) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.4.and.kx.eq.5) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.4.and.kx.eq.6) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.4.and.kx.eq.7) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.4.and.kx.eq.8) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.5.and.kx.eq.1) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.5.and.kx.eq.2) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.5.and.kx.eq.3) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.5.and.kx.eq.4) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.5.and.kx.eq.5) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.5.and.kx.eq.6) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.5.and.kx.eq.7) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.5.and.kx.eq.8) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.6.and.kx.eq.1) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.6.and.kx.eq.2) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.6.and.kx.eq.3) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.6.and.kx.eq.4) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.6.and.kx.eq.5) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.6.and.kx.eq.6) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.6.and.kx.eq.7) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.6.and.kx.eq.8) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.7.and.kx.eq.1) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.7.and.kx.eq.2) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.7.and.kx.eq.3) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.7.and.kx.eq.4) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.7.and.kx.eq.5) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.7.and.kx.eq.6) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.7.and.kx.eq.7) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.7.and.kx.eq.8) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.8.and.kx.eq.1) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.8.and.kx.eq.2) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.8.and.kx.eq.3) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.8.and.kx.eq.4) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.8.and.kx.eq.5) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.8.and.kx.eq.6) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.8.and.kx.eq.7) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.8.and.kx.eq.8) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.9.and.kx.eq.1) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.9.and.kx.eq.2) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.9.and.kx.eq.3) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.9.and.kx.eq.4) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.9.and.kx.eq.5) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.9.and.kx.eq.6) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.9.and.kx.eq.7) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.9.and.kx.eq.8) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.10.and.kx.eq.1) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.10.and.kx.eq.2) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.10.and.kx.eq.3) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.10.and.kx.eq.4) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.10.and.kx.eq.5) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.10.and.kx.eq.6) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.10.and.kx.eq.7) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.10.and.kx.eq.8) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.11.and.kx.eq.1) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.11.and.kx.eq.2) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.11.and.kx.eq.3) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.11.and.kx.eq.4) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.11.and.kx.eq.5) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.11.and.kx.eq.6) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.11.and.kx.eq.7) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.11.and.kx.eq.8) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.12.and.kx.eq.1) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.12.and.kx.eq.2) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.12.and.kx.eq.3) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.12.and.kx.eq.4) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.12.and.kx.eq.5) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.12.and.kx.eq.6) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.12.and.kx.eq.7) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.12.and.kx.eq.8) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.13.and.kx.eq.1) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.13.and.kx.eq.2) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.13.and.kx.eq.3) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.13.and.kx.eq.4) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.13.and.kx.eq.5) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.13.and.kx.eq.6) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.13.and.kx.eq.7) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.13.and.kx.eq.8) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.14.and.kx.eq.1) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.14.and.kx.eq.2) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.14.and.kx.eq.3) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.14.and.kx.eq.4) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.14.and.kx.eq.5) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.14.and.kx.eq.6) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.14.and.kx.eq.7) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.14.and.kx.eq.8) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.15.and.kx.eq.1) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.15.and.kx.eq.2) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.15.and.kx.eq.3) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.15.and.kx.eq.4) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.15.and.kx.eq.5) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.15.and.kx.eq.6) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.15.and.kx.eq.7) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.15.and.kx.eq.8) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.16.and.kx.eq.1) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.16.and.kx.eq.2) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.16.and.kx.eq.3) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.16.and.kx.eq.4) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.16.and.kx.eq.5) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.16.and.kx.eq.6) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.16.and.kx.eq.7) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.16.and.kx.eq.8) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.17.and.kx.eq.1) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.17.and.kx.eq.2) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.17.and.kx.eq.3) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.17.and.kx.eq.4) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.17.and.kx.eq.5) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.17.and.kx.eq.6) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.17.and.kx.eq.7) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.17.and.kx.eq.8) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.18.and.kx.eq.1) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.18.and.kx.eq.2) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.18.and.kx.eq.3) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.18.and.kx.eq.4) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.18.and.kx.eq.5) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.18.and.kx.eq.6) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.18.and.kx.eq.7) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.18.and.kx.eq.8) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.19.and.kx.eq.1) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.19.and.kx.eq.2) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.19.and.kx.eq.3) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.19.and.kx.eq.4) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.19.and.kx.eq.5) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.19.and.kx.eq.6) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.19.and.kx.eq.7) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.19.and.kx.eq.8) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.20.and.kx.eq.1) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.20.and.kx.eq.2) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.20.and.kx.eq.3) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.20.and.kx.eq.4) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.20.and.kx.eq.5) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.20.and.kx.eq.6) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.20.and.kx.eq.7) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                  else IF (ky.eq.20.and.kx.eq.8) then
                    ct4(kx,ky)=ct4(kx,ky)+1
                 end if                         
               end if  
               if (mevac(kx,ky,n4).eq.5) then 
                  cnt5=cnt5+1 
                  IF (ky.eq.1.and.kx.eq.1) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.1.and.kx.eq.2) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.1.and.kx.eq.3) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.1.and.kx.eq.4) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.1.and.kx.eq.5) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.1.and.kx.eq.6) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.1.and.kx.eq.7) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.1.and.kx.eq.8) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.2.and.kx.eq.1) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.2.and.kx.eq.2) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.2.and.kx.eq.3) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.2.and.kx.eq.4) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.2.and.kx.eq.5) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.2.and.kx.eq.6) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.2.and.kx.eq.7) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.2.and.kx.eq.8) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.3.and.kx.eq.1) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.3.and.kx.eq.2) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.3.and.kx.eq.3) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.3.and.kx.eq.4) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.3.and.kx.eq.5) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.3.and.kx.eq.6) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.3.and.kx.eq.7) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.3.and.kx.eq.8) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.4.and.kx.eq.1) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.4.and.kx.eq.2) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.4.and.kx.eq.3) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.4.and.kx.eq.4) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.4.and.kx.eq.5) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.4.and.kx.eq.6) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.4.and.kx.eq.7) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.4.and.kx.eq.8) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.5.and.kx.eq.1) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.5.and.kx.eq.2) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.5.and.kx.eq.3) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.5.and.kx.eq.4) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.5.and.kx.eq.5) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.5.and.kx.eq.6) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.5.and.kx.eq.7) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.5.and.kx.eq.8) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.6.and.kx.eq.1) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.6.and.kx.eq.2) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.6.and.kx.eq.3) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.6.and.kx.eq.4) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.6.and.kx.eq.5) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.6.and.kx.eq.6) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.6.and.kx.eq.7) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.6.and.kx.eq.8) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.7.and.kx.eq.1) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.7.and.kx.eq.2) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.7.and.kx.eq.3) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.7.and.kx.eq.4) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.7.and.kx.eq.5) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.7.and.kx.eq.6) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.7.and.kx.eq.7) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.7.and.kx.eq.8) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.8.and.kx.eq.1) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.8.and.kx.eq.2) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.8.and.kx.eq.3) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.8.and.kx.eq.4) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.8.and.kx.eq.5) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.8.and.kx.eq.6) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.8.and.kx.eq.7) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.8.and.kx.eq.8) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.9.and.kx.eq.1) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.9.and.kx.eq.2) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.9.and.kx.eq.3) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.9.and.kx.eq.4) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.9.and.kx.eq.5) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.9.and.kx.eq.6) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.9.and.kx.eq.7) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.9.and.kx.eq.8) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.10.and.kx.eq.1) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.10.and.kx.eq.2) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.10.and.kx.eq.3) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.10.and.kx.eq.4) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.10.and.kx.eq.5) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.10.and.kx.eq.6) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.10.and.kx.eq.7) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.10.and.kx.eq.8) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.11.and.kx.eq.1) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.11.and.kx.eq.2) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.11.and.kx.eq.3) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.11.and.kx.eq.4) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.11.and.kx.eq.5) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.11.and.kx.eq.6) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.11.and.kx.eq.7) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.11.and.kx.eq.8) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.12.and.kx.eq.1) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.12.and.kx.eq.2) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.12.and.kx.eq.3) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.12.and.kx.eq.4) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.12.and.kx.eq.5) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.12.and.kx.eq.6) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.12.and.kx.eq.7) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.12.and.kx.eq.8) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.13.and.kx.eq.1) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.13.and.kx.eq.2) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.13.and.kx.eq.3) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.13.and.kx.eq.4) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.13.and.kx.eq.5) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.13.and.kx.eq.6) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.13.and.kx.eq.7) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.13.and.kx.eq.8) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.14.and.kx.eq.1) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.14.and.kx.eq.2) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.14.and.kx.eq.3) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.14.and.kx.eq.4) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.14.and.kx.eq.5) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.14.and.kx.eq.6) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.14.and.kx.eq.7) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.14.and.kx.eq.8) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.15.and.kx.eq.1) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.15.and.kx.eq.2) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.15.and.kx.eq.3) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.15.and.kx.eq.4) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.15.and.kx.eq.5) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.15.and.kx.eq.6) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.15.and.kx.eq.7) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.15.and.kx.eq.8) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.16.and.kx.eq.1) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.16.and.kx.eq.2) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.16.and.kx.eq.3) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.16.and.kx.eq.4) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.16.and.kx.eq.5) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.16.and.kx.eq.6) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.16.and.kx.eq.7) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.16.and.kx.eq.8) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.17.and.kx.eq.1) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.17.and.kx.eq.2) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.17.and.kx.eq.3) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.17.and.kx.eq.4) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.17.and.kx.eq.5) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.17.and.kx.eq.6) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.17.and.kx.eq.7) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.17.and.kx.eq.8) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.18.and.kx.eq.1) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.18.and.kx.eq.2) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.18.and.kx.eq.3) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.18.and.kx.eq.4) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.18.and.kx.eq.5) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.18.and.kx.eq.6) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.18.and.kx.eq.7) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.18.and.kx.eq.8) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.19.and.kx.eq.1) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.19.and.kx.eq.2) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.19.and.kx.eq.3) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.19.and.kx.eq.4) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.19.and.kx.eq.5) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.19.and.kx.eq.6) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.19.and.kx.eq.7) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.19.and.kx.eq.8) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.20.and.kx.eq.1) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.20.and.kx.eq.2) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.20.and.kx.eq.3) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.20.and.kx.eq.4) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.20.and.kx.eq.5) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.20.and.kx.eq.6) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.20.and.kx.eq.7) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  else IF (ky.eq.20.and.kx.eq.8) then
                    ct5(kx,ky)=ct5(kx,ky)+1
                  end if                           
               end if 
c.... 8s are mevac=1 who are temporarily waiting for open spot.... 
c.... but not long enough to be mevac=3.....
               if (mevac(kx,ky,n4).eq.8) then 
                  cnt8=cnt8+1 
                  IF (ky.eq.1.and.kx.eq.1) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.1.and.kx.eq.2) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.1.and.kx.eq.3) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.1.and.kx.eq.4) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.1.and.kx.eq.5) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.1.and.kx.eq.6) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.1.and.kx.eq.7) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.1.and.kx.eq.8) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.2.and.kx.eq.1) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.2.and.kx.eq.2) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.2.and.kx.eq.3) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.2.and.kx.eq.4) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.2.and.kx.eq.5) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.2.and.kx.eq.6) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.2.and.kx.eq.7) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.2.and.kx.eq.8) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.3.and.kx.eq.1) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.3.and.kx.eq.2) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.3.and.kx.eq.3) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.3.and.kx.eq.4) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.3.and.kx.eq.5) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.3.and.kx.eq.6) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.3.and.kx.eq.7) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.3.and.kx.eq.8) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.4.and.kx.eq.1) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.4.and.kx.eq.2) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.4.and.kx.eq.3) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.4.and.kx.eq.4) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.4.and.kx.eq.5) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.4.and.kx.eq.6) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.4.and.kx.eq.7) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.4.and.kx.eq.8) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.5.and.kx.eq.1) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.5.and.kx.eq.2) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.5.and.kx.eq.3) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.5.and.kx.eq.4) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.5.and.kx.eq.5) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.5.and.kx.eq.6) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.5.and.kx.eq.7) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.5.and.kx.eq.8) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.6.and.kx.eq.1) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.6.and.kx.eq.2) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.6.and.kx.eq.3) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.6.and.kx.eq.4) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.6.and.kx.eq.5) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.6.and.kx.eq.6) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.6.and.kx.eq.7) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.6.and.kx.eq.8) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.7.and.kx.eq.1) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.7.and.kx.eq.2) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.7.and.kx.eq.3) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.7.and.kx.eq.4) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.7.and.kx.eq.5) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.7.and.kx.eq.6) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.7.and.kx.eq.7) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.7.and.kx.eq.8) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.8.and.kx.eq.1) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.8.and.kx.eq.2) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.8.and.kx.eq.3) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.8.and.kx.eq.4) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.8.and.kx.eq.5) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.8.and.kx.eq.6) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.8.and.kx.eq.7) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.8.and.kx.eq.8) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.9.and.kx.eq.1) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.9.and.kx.eq.2) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.9.and.kx.eq.3) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.9.and.kx.eq.4) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.9.and.kx.eq.5) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.9.and.kx.eq.6) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.9.and.kx.eq.7) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.9.and.kx.eq.8) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.10.and.kx.eq.1) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.10.and.kx.eq.2) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.10.and.kx.eq.3) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.10.and.kx.eq.4) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.10.and.kx.eq.5) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.10.and.kx.eq.6) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.10.and.kx.eq.7) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.10.and.kx.eq.8) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.11.and.kx.eq.1) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.11.and.kx.eq.2) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.11.and.kx.eq.3) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.11.and.kx.eq.4) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.11.and.kx.eq.5) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.11.and.kx.eq.6) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.11.and.kx.eq.7) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.11.and.kx.eq.8) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.12.and.kx.eq.1) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.12.and.kx.eq.2) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.12.and.kx.eq.3) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.12.and.kx.eq.4) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.12.and.kx.eq.5) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.12.and.kx.eq.6) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.12.and.kx.eq.7) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.12.and.kx.eq.8) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.13.and.kx.eq.1) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.13.and.kx.eq.2) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.13.and.kx.eq.3) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.13.and.kx.eq.4) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.13.and.kx.eq.5) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.13.and.kx.eq.6) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.13.and.kx.eq.7) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.13.and.kx.eq.8) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.14.and.kx.eq.1) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.14.and.kx.eq.2) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.14.and.kx.eq.3) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.14.and.kx.eq.4) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.14.and.kx.eq.5) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.14.and.kx.eq.6) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.14.and.kx.eq.7) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.14.and.kx.eq.8) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.15.and.kx.eq.1) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.15.and.kx.eq.2) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.15.and.kx.eq.3) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.15.and.kx.eq.4) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.15.and.kx.eq.5) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.15.and.kx.eq.6) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.15.and.kx.eq.7) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.15.and.kx.eq.8) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.16.and.kx.eq.1) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.16.and.kx.eq.2) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.16.and.kx.eq.3) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.16.and.kx.eq.4) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.16.and.kx.eq.5) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.16.and.kx.eq.6) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.16.and.kx.eq.7) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.16.and.kx.eq.8) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.17.and.kx.eq.1) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.17.and.kx.eq.2) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.17.and.kx.eq.3) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.17.and.kx.eq.4) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.17.and.kx.eq.5) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.17.and.kx.eq.6) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.17.and.kx.eq.7) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.17.and.kx.eq.8) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.18.and.kx.eq.1) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.18.and.kx.eq.2) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.18.and.kx.eq.3) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.18.and.kx.eq.4) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.18.and.kx.eq.5) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.18.and.kx.eq.6) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.18.and.kx.eq.7) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.18.and.kx.eq.8) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.19.and.kx.eq.1) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.19.and.kx.eq.2) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.19.and.kx.eq.3) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.19.and.kx.eq.4) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.19.and.kx.eq.5) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.19.and.kx.eq.6) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.19.and.kx.eq.7) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.19.and.kx.eq.8) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.20.and.kx.eq.1) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.20.and.kx.eq.2) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.20.and.kx.eq.3) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.20.and.kx.eq.4) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.20.and.kx.eq.5) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.20.and.kx.eq.6) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.20.and.kx.eq.7) then
                    ct8(kx,ky)=ct8(kx,ky)+1
                  else IF (ky.eq.20.and.kx.eq.8) then
                    ct8(kx,ky)=ct8(kx,ky)+1					
                 end if                           
               end if    
c.... 8s are mevac=1 who are temporarily waiting for open spot.... 
c.... but not long enough to be mevac=3.....
               if (mevac(kx,ky,n4).eq.99) then 
                  cnt99=cnt99+1 
                  IF (ky.eq.1.and.kx.eq.1) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.1.and.kx.eq.2) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.1.and.kx.eq.3) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.1.and.kx.eq.4) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.1.and.kx.eq.5) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.1.and.kx.eq.6) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.1.and.kx.eq.7) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.1.and.kx.eq.8) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.2.and.kx.eq.1) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.2.and.kx.eq.2) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.2.and.kx.eq.3) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.2.and.kx.eq.4) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.2.and.kx.eq.5) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.2.and.kx.eq.6) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.2.and.kx.eq.7) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.2.and.kx.eq.8) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.3.and.kx.eq.1) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.3.and.kx.eq.2) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.3.and.kx.eq.3) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.3.and.kx.eq.4) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.3.and.kx.eq.5) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.3.and.kx.eq.6) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.3.and.kx.eq.7) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.3.and.kx.eq.8) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.4.and.kx.eq.1) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.4.and.kx.eq.2) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.4.and.kx.eq.3) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.4.and.kx.eq.4) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.4.and.kx.eq.5) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.4.and.kx.eq.6) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.4.and.kx.eq.7) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.4.and.kx.eq.8) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.5.and.kx.eq.1) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.5.and.kx.eq.2) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.5.and.kx.eq.3) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.5.and.kx.eq.4) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.5.and.kx.eq.5) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.5.and.kx.eq.6) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.5.and.kx.eq.7) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.5.and.kx.eq.8) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.6.and.kx.eq.1) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.6.and.kx.eq.2) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.6.and.kx.eq.3) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.6.and.kx.eq.4) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.6.and.kx.eq.5) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.6.and.kx.eq.6) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.6.and.kx.eq.7) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.6.and.kx.eq.8) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.7.and.kx.eq.1) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.7.and.kx.eq.2) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.7.and.kx.eq.3) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.7.and.kx.eq.4) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.7.and.kx.eq.5) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.7.and.kx.eq.6) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.7.and.kx.eq.7) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.7.and.kx.eq.8) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.8.and.kx.eq.1) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.8.and.kx.eq.2) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.8.and.kx.eq.3) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.8.and.kx.eq.4) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.8.and.kx.eq.5) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.8.and.kx.eq.6) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.8.and.kx.eq.7) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.8.and.kx.eq.8) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.9.and.kx.eq.1) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.9.and.kx.eq.2) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.9.and.kx.eq.3) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.9.and.kx.eq.4) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.9.and.kx.eq.5) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.9.and.kx.eq.6) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.9.and.kx.eq.7) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.9.and.kx.eq.8) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.10.and.kx.eq.1) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.10.and.kx.eq.2) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.10.and.kx.eq.3) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.10.and.kx.eq.4) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.10.and.kx.eq.5) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.10.and.kx.eq.6) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.10.and.kx.eq.7) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.10.and.kx.eq.8) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.11.and.kx.eq.1) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.11.and.kx.eq.2) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.11.and.kx.eq.3) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.11.and.kx.eq.4) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.11.and.kx.eq.5) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.11.and.kx.eq.6) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.11.and.kx.eq.7) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.11.and.kx.eq.8) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.12.and.kx.eq.1) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.12.and.kx.eq.2) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.12.and.kx.eq.3) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.12.and.kx.eq.4) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.12.and.kx.eq.5) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.12.and.kx.eq.6) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.12.and.kx.eq.7) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.12.and.kx.eq.8) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.13.and.kx.eq.1) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.13.and.kx.eq.2) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.13.and.kx.eq.3) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.13.and.kx.eq.4) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.13.and.kx.eq.5) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.13.and.kx.eq.6) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.13.and.kx.eq.7) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.13.and.kx.eq.8) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.14.and.kx.eq.1) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.14.and.kx.eq.2) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.14.and.kx.eq.3) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.14.and.kx.eq.4) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.14.and.kx.eq.5) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.14.and.kx.eq.6) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.14.and.kx.eq.7) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.14.and.kx.eq.8) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.15.and.kx.eq.1) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.15.and.kx.eq.2) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.15.and.kx.eq.3) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.15.and.kx.eq.4) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.15.and.kx.eq.5) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.15.and.kx.eq.6) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.15.and.kx.eq.7) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.15.and.kx.eq.8) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.16.and.kx.eq.1) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.16.and.kx.eq.2) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.16.and.kx.eq.3) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.16.and.kx.eq.4) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.16.and.kx.eq.5) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.16.and.kx.eq.6) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.16.and.kx.eq.7) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.16.and.kx.eq.8) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.17.and.kx.eq.1) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.17.and.kx.eq.2) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.17.and.kx.eq.3) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.17.and.kx.eq.4) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.17.and.kx.eq.5) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.17.and.kx.eq.6) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.17.and.kx.eq.7) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.17.and.kx.eq.8) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.18.and.kx.eq.1) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.18.and.kx.eq.2) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.18.and.kx.eq.3) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.18.and.kx.eq.4) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.18.and.kx.eq.5) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.18.and.kx.eq.6) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.18.and.kx.eq.7) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.18.and.kx.eq.8) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.19.and.kx.eq.1) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.19.and.kx.eq.2) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.19.and.kx.eq.3) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.19.and.kx.eq.4) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.19.and.kx.eq.5) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.19.and.kx.eq.6) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.19.and.kx.eq.7) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.19.and.kx.eq.8) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.20.and.kx.eq.1) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.20.and.kx.eq.2) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.20.and.kx.eq.3) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.20.and.kx.eq.4) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.20.and.kx.eq.5) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.20.and.kx.eq.6) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.20.and.kx.eq.7) then
                    ct99(kx,ky)=ct99(kx,ky)+1
                  else IF (ky.eq.20.and.kx.eq.8) then
                    ct99(kx,ky)=ct99(kx,ky)+1                   
                 end if                           
               end if                   
               ct6(kx,ky)=ct0(kx,ky)+ct3(kx,ky)+ct4(kx,ky)
     *            +ct5(kx,ky)+ct99(kx,ky)
               ct9(kx,ky)=ct1(kx,ky)+ct8(kx,ky)
               ct10(kx,ky)=ct6(kx,ky)+ct9(kx,ky)+ct2(kx,ky)
               ct11(kx,ky)=ct1(kx,ky)+ct2(kx,ky)+ct8(kx,ky)
               ct12(kx,ky)=ct11(kx,ky)+ct6(kx,ky)
             END DO
c             end if
           end do
         end do
 333     continue
          ncount=ncount+1
          if (ncount.eq.50) then
            ncount=0
            nmin=nmin+1
            if (mod(nmin,5).eq.0) write(*,*) 'completed min: ',nmin
          end if
c. This section will PRINT out the TRAFFIC into various large N Files
c... print out between 2h-2.1h..... EW Interstates 
c. its commented out currently but this is how its done
c          if (ntick.gt.16000.and.ntick.le.16003) then
c            do i=1,NCELLY
c            do ln=1,7
c              ky=1
c                write(69,'(a1,1x,i6,1x,i6,1x,i6,1x,i5,1x,i8,4(1x,i1))') 
c     *             'E',ky,ntick,i,ln,ident(i,ln),
c*             car(i,ln),vel(i,ln),twoe(i,ln),exite(i,ln)
c                write(69,'(a1,1x,i6,1x,i6,1x,i6,1x,i5,1x,i8,4(1x,i1))') 
c     *             'W',ky,ntick,i,ln,identX(i,ln),
c     *             carX(i,ln),velX(i,ln),twow(i,ln),exitw(i,ln)
c           end do
c            end do
c          end if
c          
c          if (ntick.gt.1000.and.ntick.le.1200) then
c            ky=10
c            kx=3
c            do i=1,NCELLX
c            do lx=1,2
c                write(21,'(a1,1x,i6,1x,i6,1x,i6,1x,i5,1x,i8,2(1x,i1))') 
c     *              '1',ky,ntick,i,lx,idente(kx,ky,i,lx),
c     *                care(kx,ky,i,lx),vele(kx,ky,i,lx)
c            end do   
c            END DO
c          end if     

c. here we print out - both to the screen and different fort files...
c. information on the number of evacuees etc. summed up across grid cells 
c. this is done every 6 hours 
        
          if (mod(ntick,18000).eq.0) then
            write(14,*) nmin/60,cnt0*4,cnt1*4,cnt2*4,cnt3*4,cnt4*4,
     *               cnt5*4,cnt8*4,cnt99*4
            write(15,*) 'nmin',nmin,'hour',nmin/60
            write(15,*) 'spots avail, used'
            do ly=21,4,-1
              write(15,*) ly,
     *          ((maxspots(lx,ly)-spt(lx,ly))*4,lx=1,8),
     *          (spt(lx,ly)*4,lx=1,8)
            end do
            write(15,*) '%evac //%left//%wait//%gave up'            
            do ky=20,1,-1
             write(15,*) ky,(((ct11(kx,ky)*100)/ct12(kx,ky)),kx=1,8),
     *          ((ct2(kx,ky)*100)/ct10(kx,ky),kx=1,8),
     *          ((ct9(kx,ky)*100)/ct10(kx,ky),kx=1,8),
     *          ((ct3(kx,ky)*100)/ct10(kx,ky),kx=1,8)            
            end do
            write(15,*) '%notdeciding'            
            do ky=20,1,-1
             write(15,*) ky,(((ct99(kx,ky)*100)/ct10(kx,ky)),kx=1,8)           
            end do
          end if 

c. this generates some files at different time periods that i convert to csv files
c. and visualize the evacuation over time in python

           if (ntick.eq.(18000*2)) then
             do ky=1,10
                do kx=1,8              
                   write(27,*) kx,ky,((ct2(kx,ky)*100)/ct10(kx,ky)),
     *          	             ((ct3(kx,ky)*100)/ct10(kx,ky))				
                end do
             end do 
             do ky=11,20
                do kx=1,8              
                   write(27,*) kx,ky,((ct2(kx,ky)*100)/ct10(kx,ky)),
     *          	             ((ct3(kx,ky)*100)/ct10(kx,ky))				
                end do
             end do   
           end if
           if (ntick.eq.(18000*3)) then
             do ky=1,10
                do kx=1,8              
                   write(28,*) kx,ky,((ct2(kx,ky)*100)/ct10(kx,ky)),
     *          	             ((ct3(kx,ky)*100)/ct10(kx,ky))				
                end do
             end do 
             do ky=11,20
                do kx=1,8              
                   write(28,*) kx,ky,((ct2(kx,ky)*100)/ct10(kx,ky)),
     *          	             ((ct3(kx,ky)*100)/ct10(kx,ky))				
                end do
             end do    
           end if
           if (ntick.eq.(18000*4)) then
             do ky=1,10
                do kx=1,8              
                   write(29,*) kx,ky,((ct2(kx,ky)*100)/ct10(kx,ky)),
     *          	             ((ct3(kx,ky)*100)/ct10(kx,ky))				
                end do
             end do 
             do ky=11,20
                do kx=1,8              
                   write(29,*) kx,ky,((ct2(kx,ky)*100)/ct10(kx,ky)),
     *          	             ((ct3(kx,ky)*100)/ct10(kx,ky))				
                end do
             end do   
           end if
           if (ntick.eq.(18000*5)) then
             do ky=1,10
                do kx=1,8              
                   write(30,*) kx,ky,((ct2(kx,ky)*100)/ct10(kx,ky)),
     *          	             ((ct3(kx,ky)*100)/ct10(kx,ky))				
                end do
             end do 
             do ky=11,20
                do kx=1,8              
                   write(30,*) kx,ky,((ct2(kx,ky)*100)/ct10(kx,ky)),
     *          	             ((ct3(kx,ky)*100)/ct10(kx,ky))				
                end do
             end do   
           end if
           if (ntick.eq.(18000*6)) then
             do ky=1,10
                do kx=1,8              
                   write(31,*) kx,ky,((ct2(kx,ky)*100)/ct10(kx,ky)),
     *          	             ((ct3(kx,ky)*100)/ct10(kx,ky))				
                end do
             end do 
             do ky=11,20
                do kx=1,8              
                   write(31,*) kx,ky,((ct2(kx,ky)*100)/ct10(kx,ky)),
     *          	             ((ct3(kx,ky)*100)/ct10(kx,ky))				
                end do
             end do   
           end if
		     if (ntick.eq.(18000*7)) then
             do ky=1,10
                do kx=1,8              
                   write(32,*) kx,ky,((ct2(kx,ky)*100)/ct10(kx,ky)),
     *          	             ((ct3(kx,ky)*100)/ct10(kx,ky))				
                end do
             end do 
             do ky=11,20
                do kx=1,8              
                   write(32,*) kx,ky,((ct2(kx,ky)*100)/ct10(kx,ky)),
     *          	             ((ct3(kx,ky)*100)/ct10(kx,ky))				
                end do
             end do   
           end if
           if (ntick.eq.(18000*8)) then
             do ky=1,10
                do kx=1,8              
                   write(33,*) kx,ky,((ct2(kx,ky)*100)/ct10(kx,ky)),
     *          	             ((ct3(kx,ky)*100)/ct10(kx,ky))				
                end do
             end do 
             do ky=11,20
                do kx=1,8              
                   write(33,*) kx,ky,((ct2(kx,ky)*100)/ct10(kx,ky)),
     *          	             ((ct3(kx,ky)*100)/ct10(kx,ky))				
                end do
             end do   
           end if
           if (ntick.eq.(18000*9)) then
             do ky=1,10
                do kx=1,8              
                   write(34,*) kx,ky,((ct2(kx,ky)*100)/ct10(kx,ky)),
     *          	             ((ct3(kx,ky)*100)/ct10(kx,ky))				
                end do
             end do 
             do ky=11,20
                do kx=1,8              
                   write(34,*) kx,ky,((ct2(kx,ky)*100)/ct10(kx,ky)),
     *          	             ((ct3(kx,ky)*100)/ct10(kx,ky))				
                end do
             end do   
           end if
           if (ntick.eq.(18000*10)) then
             do ky=1,10
                do kx=1,8              
                   write(35,*) kx,ky,((ct2(kx,ky)*100)/ct10(kx,ky)),
     *          	             ((ct3(kx,ky)*100)/ct10(kx,ky))				
                end do
             end do 
             do ky=11,20
                do kx=1,8              
                   write(35,*) kx,ky,((ct2(kx,ky)*100)/ct10(kx,ky)),
     *          	             ((ct3(kx,ky)*100)/ct10(kx,ky))				
                end do
             end do   
           end if
           if (ntick.eq.(18000*11)) then
             do ky=1,10
                do kx=1,8              
                   write(36,*) kx,ky,((ct2(kx,ky)*100)/ct10(kx,ky)),
     *          	             ((ct3(kx,ky)*100)/ct10(kx,ky))				
                end do
             end do 
             do ky=11,20
                do kx=1,8              
                   write(36,*) kx,ky,((ct2(kx,ky)*100)/ct10(kx,ky)),
     *          	             ((ct3(kx,ky)*100)/ct10(kx,ky))				
                end do
             end do    
           end if
           if (ntick.eq.(18000*12)) then
             do ky=1,10
                do kx=1,8              
                   write(37,*) kx,ky,((ct2(kx,ky)*100)/ct10(kx,ky)),
     *          	             ((ct3(kx,ky)*100)/ct10(kx,ky))				
                end do
             end do 
             do ky=11,20
                do kx=1,8              
                   write(37,*) kx,ky,((ct2(kx,ky)*100)/ct10(kx,ky)),
     *          	             ((ct3(kx,ky)*100)/ct10(kx,ky))				
                end do
             end do   
           end if
           if (ntick.eq.(18000*13)) then
             do ky=1,10
                do kx=1,8              
                   write(38,*) kx,ky,((ct2(kx,ky)*100)/ct10(kx,ky)),
     *          	             ((ct3(kx,ky)*100)/ct10(kx,ky))				
                end do
             end do 
             do ky=11,20
                do kx=1,8              
                   write(38,*) kx,ky,((ct2(kx,ky)*100)/ct10(kx,ky)),
     *          	             ((ct3(kx,ky)*100)/ct10(kx,ky))				
                end do
             end do   
           end if
           if (ntick.eq.(18000*14)) then
             do ky=1,10
                do kx=1,8              
                   write(39,*) kx,ky,((ct2(kx,ky)*100)/ct10(kx,ky)),
     *          	             ((ct3(kx,ky)*100)/ct10(kx,ky))				
                end do
             end do 
             do ky=11,20
                do kx=1,8              
                   write(39,*) kx,ky,((ct2(kx,ky)*100)/ct10(kx,ky)),
     *          	             ((ct3(kx,ky)*100)/ct10(kx,ky))				
                end do
             end do   
           end if
           if (ntick.eq.(18000*15)) then
             do ky=1,10
                do kx=1,8              
                   write(40,*) kx,ky,((ct2(kx,ky)*100)/ct10(kx,ky)),
     *          	             ((ct3(kx,ky)*100)/ct10(kx,ky))				
                end do
             end do 
             do ky=11,20
                do kx=1,8              
                   write(40,*) kx,ky,((ct2(kx,ky)*100)/ct10(kx,ky)),
     *          	             ((ct3(kx,ky)*100)/ct10(kx,ky))				
                end do
             end do   
           end if
           if (ntick.eq.(18000*16)) then
             do ky=1,10
                do kx=1,8              
                   write(41,*) kx,ky,((ct2(kx,ky)*100)/ct10(kx,ky)),
     *          	             ((ct3(kx,ky)*100)/ct10(kx,ky))				
                end do
             end do 
             do ky=11,20
                do kx=1,8              
                   write(41,*) kx,ky,((ct2(kx,ky)*100)/ct10(kx,ky)),
     *          	             ((ct3(kx,ky)*100)/ct10(kx,ky))				
                end do
             end do   
           end if
           if (ntick.eq.(18000*17)) then
             do ky=1,10
                do kx=1,8              
                   write(42,*) kx,ky,((ct2(kx,ky)*100)/ct10(kx,ky)),
     *          	             ((ct3(kx,ky)*100)/ct10(kx,ky))				
                end do
             end do 
             do ky=11,20
                do kx=1,8              
                   write(42,*) kx,ky,((ct2(kx,ky)*100)/ct10(kx,ky)),
     *          	             ((ct3(kx,ky)*100)/ct10(kx,ky))				
                end do
             end do   
           end if
           if (ntick.eq.(18000*18)) then
             do ky=1,10
                do kx=1,8              
                   write(43,*) kx,ky,((ct2(kx,ky)*100)/ct10(kx,ky)),
     *          	             ((ct3(kx,ky)*100)/ct10(kx,ky))				
                end do
             end do 
             do ky=11,20
                do kx=1,8              
                   write(43,*) kx,ky,((ct2(kx,ky)*100)/ct10(kx,ky)),
     *          	             ((ct3(kx,ky)*100)/ct10(kx,ky))				
                end do
             end do   
           end if
           if (ntick.eq.(18000*19)) then
             do ky=1,10
                do kx=1,8              
                   write(44,*) kx,ky,((ct2(kx,ky)*100)/ct10(kx,ky)),
     *          	             ((ct3(kx,ky)*100)/ct10(kx,ky))				
                end do
             end do 
             do ky=11,20
                do kx=1,8              
                   write(44,*) kx,ky,((ct2(kx,ky)*100)/ct10(kx,ky)),
     *          	             ((ct3(kx,ky)*100)/ct10(kx,ky))				
                end do
             end do   
           end if
           if (ntick.eq.(18000*20)) then
             do ky=1,10
                do kx=1,8              
                   write(45,*) kx,ky,((ct2(kx,ky)*100)/ct10(kx,ky)),
     *          	             ((ct3(kx,ky)*100)/ct10(kx,ky))				
                end do
             end do 
             do ky=11,20
                do kx=1,8              
                   write(45,*) kx,ky,((ct2(kx,ky)*100)/ct10(kx,ky)),
     *          	             ((ct3(kx,ky)*100)/ct10(kx,ky))				
                end do
             end do    
           end if
           if (ntick.eq.(18000*21)) then
             do ky=1,10
                do kx=1,8              
                   write(46,*) kx,ky,((ct2(kx,ky)*100)/ct10(kx,ky)),
     *          	             ((ct3(kx,ky)*100)/ct10(kx,ky))				
                end do
             end do 
             do ky=11,20
                do kx=1,8              
                   write(46,*) kx,ky,((ct2(kx,ky)*100)/ct10(kx,ky)),
     *          	             ((ct3(kx,ky)*100)/ct10(kx,ky))				
                end do
             end do   
           end if
           if (ntick.eq.(18000*22)) then
             do ky=1,10
                do kx=1,8              
                   write(47,*) kx,ky,((ct2(kx,ky)*100)/ct10(kx,ky)),
     *          	             ((ct3(kx,ky)*100)/ct10(kx,ky))				
                end do
             end do 
             do ky=11,20
                do kx=1,8              
                   write(47,*) kx,ky,((ct2(kx,ky)*100)/ct10(kx,ky)),
     *          	             ((ct3(kx,ky)*100)/ct10(kx,ky))				
                end do
             end do   
           end if
           if (ntick.eq.(18000*23)) then
             do ky=1,10
                do kx=1,8              
                   write(48,*) kx,ky,((ct2(kx,ky)*100)/ct10(kx,ky)),
     *          	             ((ct3(kx,ky)*100)/ct10(kx,ky))				
                end do
             end do 
             do ky=11,20
                do kx=1,8              
                   write(48,*) kx,ky,((ct2(kx,ky)*100)/ct10(kx,ky)),
     *          	             ((ct3(kx,ky)*100)/ct10(kx,ky))				
                end do
             end do    
           end if
		   
c. This section will print out an update to the screen every hour
c. these are total stats across all grid cells
          if (mod(ntick,3000).eq.0) then
            write(*,*) 'Total pop in risk zones'            
            write(*,*) (cnt0+cnt1+cnt2+cnt3+cnt4+cnt5+cnt8+cnt99)*4
            write(*,*) 'pop staying (mevac=0)*4'            
            write(*,*) cnt0*4
            write(*,*) 'pop wanting to leave (mevac=1)*4'            
            write(*,*) cnt1*4    
            write(*,*) 'Groups left (mevac=2)*4'            
            write(*,*) cnt2*4    
            write(*,*) 'Groups who gave up, now staying (mevac=3)*4'            
            write(*,*) cnt3*4  
            write(*,*) 'Cant leave (mevac=4)*4'            
            write(*,*) cnt4*4 
            write(*,*) 'Local Shelter (mevac=5)*4'            
            write(*,*) cnt5*4     
            write(*,*) 'Looking for spot (mevac=8)'            
            write(*,*) cnt8*4  
            write(*,*) 'Not deciding (mevac=99)'            
            write(*,*) cnt99*4  
            write(10,*) 'Total pop in risk zones'            
            write(10,*) (cnt0+cnt1+cnt2+cnt3+cnt4+cnt5+cnt8+cnt99)*4
            write(10,*) 'pop staying (mevac=0)*4'            
            write(10,*) cnt0*4
            write(10,*) 'pop wanting to leave (mevac=1)*4'            
            write(10,*) cnt1*4    
            write(10,*) 'Groups left (mevac=2)*4'            
            write(10,*) cnt2*4    
            write(10,*) 'Groups who gave up, now staying (mevac=3)*4'            
            write(10,*) cnt3*4  
            write(10,*) 'Cant leave (mevac=4)*4'            
            write(10,*) cnt4*4 
            write(10,*) 'Local Shelter (mevac=5)*4'            
            write(10,*) cnt5*4     
            write(10,*) 'Looking for spot (mevac=8)'            
            write(10,*) cnt8*4    
            write(10,*) 'Not deciding (mevac=99)'            
            write(10,*) cnt99*4              
            write(*,*) ntick,nevac,nmin
            write(10,*) ntick,nevac,nmin
            write(*,*) 'completed min: ',nmin
            write(10,*) 'completed min: ',nmin 
            write(*,*) 'out of state', state*4
            write(*,*) 'max spots, spots available, spots taken,%filled'
            write(10,*) 'max spots, spots avail, spots taken,%filled'
            write(12,*) nmin,num,
     *         int((cnt0+cnt1+cnt2+cnt3+cnt4+cnt5+cnt8+cnt99)*4),
     *         int(cnt0*4),int(cnt1*4),int(cnt2*4),int(cnt3*4),
     *         int(cnt4*4),int(cnt5*4),int(cnt8*4),int(cnt99*4)
            do ly=21,8,-1
              write(*,*) ly,(maxspots(lx,ly)*4,lx=1,8),
     *          ((maxspots(lx,ly)-spt(lx,ly))*4,lx=1,8)
              write(10,*) ly,(maxspots(lx,ly)*4,lx=1,8),
     *          ((maxspots(lx,ly)-spt(lx,ly))*4,lx=1,8)
            end do
           write(*,*) 'spots taken'
            write(10,*) 'spots taken'
            do ly=21,8,-1
              write(*,*) ly,
     *          (spt(lx,ly)*4,lx=1,8)
              write(10,*) ly,
     *          (spt(lx,ly)*4,lx=1,8)
            end do
            write(*,*) 'Staying // left'
            write(10,*) 'Staying // left'
            nsum=0
            do ky=20,1,-1
              write(*,*) ky,(ct6(kx,ky)*4,kx=1,8),
     *          (ct2(kx,ky)*4,kx=1,8)
              write(10,*) ky,(ct6(kx,ky)*4,kx=1,4),
     *          (ct2(kx,ky)*4,kx=1,8)
c ... figure left at risk
              if (nmin.eq.11520.and.ivac.eq.1) then
                if (ky.le.20) then
                  do kx=1,8
                    nsum=nsum+npop(kx,ky)
                  end do
                end if
              end if
            end do
            write(*,*) 'look4spot // gave up'
            write(10,*) 'look4spot // gave up'
            nsum=0
            do ky=20,1,-1
              write(*,*) ky,
     *          (ct8(kx,ky)*4,kx=1,8), 
     *          (ct3(kx,ky)*4,kx=1,8)
              write(10,*) ky,
     *          (ct8(kx,ky)*4,kx=1,8),
     *          (ct3(kx,ky)*4,kx=1,8)
c ... figure left at risk
              if (nmin.eq.11520.and.ivac.eq.1) then
                if (ky.le.20) then
                  do kx=1,8
                    nsum=nsum+npop(kx,ky)
                  end do
                end if
              end if
            end do
            nsum=0
              write(*,*) '%evac // %left'
              write(10,*) '%evac // %left'            
            do ky=20,1,-1
              write(*,*) ky,(((ct11(kx,ky)*100)/ct12(kx,ky)),kx=1,8),
     *          ((ct2(kx,ky)*100)/ct10(kx,ky),kx=1,8)
            write(10,*) ky,(((ct11(kx,ky)*100)/ct12(kx,ky)),kx=1,8),
     *          ((ct2(kx,ky)*100)/ct10(kx,ky),kx=1,8)           
            end do
              write(*,*) '%wait // %gave up'
              write(10,*) '%wait // %gave up'            
            do ky=20,1,-1
              write(*,*) ky,
     *          ((ct9(kx,ky)*100)/ct10(kx,ky),kx=1,8), 
     *          ((ct3(kx,ky)*100)/ct10(kx,ky),kx=1,8)
            write(10,*) ky,
     *          ((ct9(kx,ky)*100)/ct10(kx,ky),kx=1,8),
     *          ((ct3(kx,ky)*100)/ct10(kx,ky),kx=1,8)            
            end do
              write(*,*) '%notdecidingyet, %decidetostay'
              write(10,*) '%notdecidingyet,%decidetostay'            
            do ky=20,1,-1
              write(*,*) ky,(((ct99(kx,ky)*100)/ct10(kx,ky)),kx=1,8),
     *          ((ct0(kx,ky)*100)/ct10(kx,ky),kx=1,8)
              write(10,*) ky,(((ct99(kx,ky)*100)/ct10(kx,ky)),kx=1,8),     
     *          ((ct0(kx,ky)*100)/ct10(kx,ky),kx=1,8)              
            end do            
c ... still on road
            numE=0
            numW=0
            numSW=0
			numSE=0
			numMW=0
			numH=0
			
            do ln=1,5
            do i=1,NCELLY
              numE=numE+carp(i,ln)
            end do
            end do
			write(*,*) 'Still enroute on I95: ',numE
            write(10,*) 'Still enrouteon I95: ',numE
			write(*,*) '% congestions on I95: ',(numE/NCELLY)
            write(10,*) '% congestions on I95: ',(numE/NCELLY)			
	 
            do ln=1,5
            do i=1,NCELLY
              numW=numW+carpX(i,ln)
            end do
            end do
			write(*,*) 'Still enroute on I75: ',numW
            write(10,*) 'Still enrouteon I75: ',numW
			write(*,*) '% congestions on I75: ',(numW/NCELLY)
            write(10,*) '% congestions on I75: ',(numW/NCELLY)

            write(19,*) nmin/60,numE,numW
	 
            lx=1
            do kx=1,NX
            do ky=2,NY
              do i=1,NCELLX
                numH=numH+carpe(kx,ky,i,lx)
              end do
            end do
            end do
			write(*,*) 'Still enroute on highways: ',numH
            write(10,*) 'Still enrouteon highways: ',numH
            
            do le=1,3
            ky=1
            do kx=1,NX
              do i=1,NCELLX
                numSE=numSE+carpse(kx,ky,i,le)
              end do
            end do            
            end do 
			write(*,*) 'Still enroute on SE: ',numSE
            write(10,*) 'Still enrouteon SE: ',numSE

            do le=1,3
            ky=1
            do kx=1,NX
              do i=1,NCELLX
                numSW=numSW+carpsw(kx,ky,i,le)
              end do
            end do            
            end do 
			write(*,*) 'Still enroute on SW: ',numSW
            write(10,*) 'Still enrouteon SW: ',numSW
			
			
            DO le=1,3
            ky=10
            do kx=1,8
              do i=1,NCELLX
                numMW=numMW+carpmw(kx,ky,i,le)
              end do
            end do   
			END DO 
            write(*,*) 'Still enroute on MW: ',numMW
            write(10,*) 'Still enroute on MW: ',numMW
            nevac=0
          end if
          if (nmin.eq.11520) goto 102

       end do

 102   continue
       stop
       end
       subroutine ABM(jinit,npop2,wrisk,srisk,rrisk,mevac,
     *     pevac,socio,age,nocar,mobl,agewt,
     *     wwt,swt,rwt,barr,mobwt,strmwt,ewt,tstoa)  
       PARAMETER (NCELLX=1750,NCELLY=35000,NX=8,NY=20,nVMX=5,nVMX2=3,
     *   RLC=7.5,DT=1.2,RP=0.35,RA=0.00005)

c... Each grouping has the following attributes:
c   Identifier (this is unique to each group and is used for tracking/updating)
c   Links (list of group identity connections – list other identifiers)
c
c... All of the above feeds into go/no-go decision. Each decision is made every 30 minutes when the ABM is called. 
c   mevac = 0 when the decision is to stay
c   mevac = 1 when the decision is to go. 
c   mevac = 2 when the group has left
c   mevac = 3 when the group has tried to leave but gave up (wait too long)
c   mevac = 4 are those that simply cannot leave due to lack of car etc. 
c   mevac = 5 are those who evac but go to LOCAL SHELTER
c   mevac = 8 are those that are going but looking for spot on road

c ... information for evacuation decision-making algorithm 

c
       integer npop2(NX,NY)
c. mevac is a variable here that tracks the evacuation status 	   
       integer mevac(NX,NY,295355)
c. these values come from the npop file 	   
       integer socio(NX,NY),age(NX,NY),nocar(NX,NY),mobl(NX,NY)  
       integer socios(NX,NY,295355),ages(NX,NY,295355)
       integer mobls(NX,NY,295355)
       integer wrisk(NX,NY),srisk(NX,NY),rrisk(NX,NY)
c. here are some variables that go into the evacuation decision algorithm	   
       integer pevac(NX,NY,295355)
       integer storm(NX,NY,295355),agewt(NX,NY,295355)
       integer wwt(NX,NY,295355),swt(NX,NY,295355),rwt(NX,NY,295355) 
       integer barr(NX,NY,295355),mobwt(NX,NY,295355)     
       integer strmwt(NX,NY,295355),ewt(NX,NY,295355)
       integer tstoa(NX,NY)
c
c. This section assigns group characteristics THE FIRST CALL ONLY
       if (jinit.eq.1) then
         do n1=1,NX
         do n2=1,NY
c ... set the characteristics for a group of 4
c ... aka each group will have one common set of char
           n3e=npop2(n1,n2)/4
           do n3=1,n3e
             mevac(n1,n2,n3)=99
c .... forecast weights  
c .. weights are informed by literature and various sensitivity analyses           
             wwt(n1,n2,n3)=1+((rand()*0.9)+0.10)*100
             swt(n1,n2,n3)=1+rand()*100
             rwt(n1,n2,n3)=1+(rand()*0.9)*100 
c... mobile home weights
              mobwt(n1,n2,n3)=1+(rand()*0.05)*100
c..... assigning weights for age
c              agewt(n1,n2,n3)=0           
              agewt(n1,n2,n3)=1+(rand()*0.01)*100             
c..... assigning weights for prisk
              strmwt(n1,n2,n3)=1+(rand()*0.8)*100
              ewt(n1,n2,n3)=1+(rand())*100    
c... assigning those who live in mobile homes versus not. percentages informed by spatial social data               
             if (mobl(n1,n2).eq.1) then 
                if (rand().le.0.05) then
                   mtemp=1
                else 
                   mtemp=0
                END IF
             end if               
             if (mobl(n1,n2).eq.2) then 
                if (rand().le.0.10) then
                   mtemp=1
                else 
                   mtemp=0
                END IF
             end if               
             if (mobl(n1,n2).eq.3) then 
                if (rand().le.0.20) then
                   mtemp=1
                else 
                   mtemp=0
                END IF
             end if  
             if (mobl(n1,n2).eq.4) then 
                if (rand().le.0.33) then
                   mtemp=1
                else 
                   mtemp=0
                END IF
             end if
             if (mobl(n1,n2).eq.5) then 
                if (rand().le.0.46) then
                   mtemp=1
                else 
                   mtemp=0
                END IF
             end if             
c.... combining mobile home tile and weights
             mobls(n1,n2,n3)=mtemp*mobwt(n1,n2,n3)
c.... creating household barriers based on social vulnerability data by grid cell             
              if (socio(n1,n2).eq.1) then
                 barr(n1,n2,n3)=((rand())+0.06)*100 
              else if (socio(n1,n2).eq.2) then
                 barr(n1,n2,n3)=((rand())+0.07)*100                   
              else if (socio(n1,n2).eq.3) then
                 barr(n1,n2,n3)=((rand())+0.08)*100                   
              else if (socio(n1,n2).eq.4) then
                 barr(n1,n2,n3)=((rand())+0.09)*100                    
              else if (socio(n1,n2).eq.5) then
                 barr(n1,n2,n3)=((rand())+0.10)*100
              end if 
c.... assigning car ownership to agents            
c.... those without a car cannot evacuate             
             if (nocar(n1,n2).eq.1) then 
                if (rand().le.0.04) then
                   mevac(n1,n2,n3)=4
                else 
                   mevac(n1,n2,n3)=99
                END IF
             end if
              if (nocar(n1,n2).eq.2) then 
                if (rand().le.0.06) then
                   mevac(n1,n2,n3)=4
                else 
                   mevac(n1,n2,n3)=99
                END IF
             end if
             if (nocar(n1,n2).eq.3) then 
                if (rand().le.0.07) then
                   mevac(n1,n2,n3)=4
                else 
                   mevac(n1,n2,n3)=99
                END IF
             end if
             if (nocar(n1,n2).eq.4) then 
                if (rand().le.0.09) then
                   mevac(n1,n2,n3)=4
                else 
                   mevac(n1,n2,n3)=99
                END IF
             end if
             if (nocar(n1,n2).eq.5) then 
                if (rand().le.0.11) then
                   mevac(n1,n2,n3)=4
                else 
                   mevac(n1,n2,n3)=99
                END IF
             end if                           
             nump=1+int(5.*rand())
           end do
         end do
         end do
c ... now construct links
       else
c....................................................................................................
c. This code is called every 30 minutes
c ... now we assign go/no go decisions 
c ... but NOT if mevac=3 or mevac = 3 or mevac = 1 already. 
c... only if people are at home and want to change their mind (mevac=0)
         do n1=1,NX
         do n2=1,NY
           n3e=npop2(n1,n2)/4 
           do n3=1,n3e  
           if (mevac(n1,n2,n3).eq.0) then          
              if (tstoa(n1,n2).ge.120) then
              barrs=barr(n1,n2,n3)+30
c. time component exists for barriers... i.e., gets smaller over time 			  
              end if       
              if (tstoa(n1,n2).lt.120.and.tstoa(n1,n2).ge.114) then
              barrs=barr(n1,n2,n3)+26
              end if                
              if (tstoa(n1,n2).lt.114.and.tstoa(n1,n2).ge.108) then
              barrs=barr(n1,n2,n3)+22
              end if    
              if (tstoa(n1,n2).lt.108.and.tstoa(n1,n2).ge.102) then
              barrs=barr(n1,n2,n3)+18
              end if    
              if (tstoa(n1,n2).lt.102.and.tstoa(n1,n2).ge.96) then
              barrs=barr(n1,n2,n3)+16
              end if    
              if (tstoa(n1,n2).lt.96.and.tstoa(n1,n2).ge.90) then
              barrs=barr(n1,n2,n3)+14
              end if    
              if (tstoa(n1,n2).lt.90.and.tstoa(n1,n2).ge.84) then
              barrs=barr(n1,n2,n3)+12
              end if    
              if (tstoa(n1,n2).lt.84.and.tstoa(n1,n2).ge.78) then
              barrs=barr(n1,n2,n3)+10
              end if    
              if (tstoa(n1,n2).lt.78.and.tstoa(n1,n2).ge.72) then
              barrs=barr(n1,n2,n3)+8
              end if    
              if (tstoa(n1,n2).lt.72.and.tstoa(n1,n2).ge.66) then
              barrs=barr(n1,n2,n3)+6
              end if    
              if (tstoa(n1,n2).lt.66.and.tstoa(n1,n2).ge.60) then
              barrs=barr(n1,n2,n3)+4
              end if                  
              if (tstoa(n1,n2).lt.60.and.tstoa(n1,n2).ge.54) then
              barrs=barr(n1,n2,n3)+2
              end if
              if (tstoa(n1,n2).lt.54.and.tstoa(n1,n2).ge.48) then
              barrs=barr(n1,n2,n3)
              end if
              if (tstoa(n1,n2).lt.48.and.tstoa(n1,n2).ge.42) then
              barrs=barr(n1,n2,n3)
              end if	
              if (tstoa(n1,n2).lt.42.and.tstoa(n1,n2).ge.36) then
              barrs=barr(n1,n2,n3)
              end if				  
              if (tstoa(n1,n2).lt.36) then
              barrs=barr(n1,n2,n3)
              end if     
             mob=mobls(n1,n2,n3)
             if (wrisk(n1,n2).gt.9) then
                  mob=mob*0.9
             end if
             if (wrisk(n1,n2).gt.8.and.wrisk(n1,n2).le.9) then
                  mob=mob*0.8
             end if
             if (wrisk(n1,n2).gt.7.and.wrisk(n1,n2).le.8) then
                  mob=mob*0.7
             end if    
             if (wrisk(n1,n2).gt.6.and.wrisk(n1,n2).le.7) then
                  mob=mob*0.6
             end if
             if (wrisk(n1,n2).gt.5.and.wrisk(n1,n2).le.6) then
                  mob=mob*0.5
             end if
             if (wrisk(n1,n2).gt.4.and.wrisk(n1,n2).le.5) then
                  mob=mob*0.4
             end if 
             if (wrisk(n1,n2).gt.3.and.wrisk(n1,n2).le.4) then
                  mob=mob*0.3
             end if
             if (wrisk(n1,n2).gt.2.and.wrisk(n1,n2).le.3) then
                  mob=mob*0.2
             end if
             if (wrisk(n1,n2).lt.2) then
                  mob=mob*0.1
             end if 	 
             prisk=0
c. wind, surge, and rain risk... combining the forecast values by grid cell with weights 			 
             w=((wrisk(n1,n2)-1)*0.4*33*(wwt(n1,n2,n3)))*0.01
             s=((srisk(n1,n2)-1)*0.4*33*(swt(n1,n2,n3)))*0.01
             r=((rrisk(n1,n2)-1)*0.4*33*(rwt(n1,n2,n3)))*0.01
             haz=max(w,s,r)
c. combining risk with warning info and household characteristics 	
c mevac=1 means evacuation decision made to leave, 0 is decision to stay 		 
             a=age(n1,n2)*20*(agewt(n1,n2,n3)*0.01)
             prisk=(haz*strmwt(n1,n2,n3)*0.01)+
     *          (pevac(n1,n2,n3)*ewt(n1,n2,n3))+
     *           mob+a
             if (mevac(n1,n2,n3).ne.4) then        
               if (prisk+int((rand()*10)).ge.
     *              barrs) then
                     mevac(n1,n2,n3)=1
               else
                  mevac(n1,n2,n3)=0               
               end if
             end if 
c            write(10,*) n1,n2,n3,mob,wrisk(n1,n2)
           END IF   
           end do
         end do
         end do
       end if
       return
       end
