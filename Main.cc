// Prog.cc
// This Program is

#include "Main.h"
int main (int par_exec1, char* pars_exec[])
{    
//Objects creation
	Mathematics * Math = new Mathematics();
    Physics * Phys = new Physics();
    
	FILE * outfile;
	outfile = fopen("Positions.dat","w");
	
	FILE * outfile2;
	outfile2 = fopen("larmor.dat","w");

    printf("*******PROPAGATION OF COSMIC RAYS*******\n\n");
    
//   v0[0]=c/3,v0[1]=c/30,v0[2]=0;
//   B0[0]=0,B0[1]=0,B0[2]=1e-4;
    m_0=m;
    
    E=1.00*pow(10,14);
    
  	Phys->Galactic(B0,r0);
 	
  	//Making the velocity with energy
  	
  	modv=Phys->Ener2vel(E,modv,c,m);
  	
//  printf("modv=%.3e \n",modv);
  	
  	Phys->SortVel(modv,m,v0); //This gives the function v0[3];
   
	printf("v0x=%.6e, v0y=%.6e, v0z=%.6e \n",v0[0],v0[1],v0[2]);

	//Modules and larmor:
	
    modv=Math->Mod(v0);	
	printf("modv0=%.10e m/s = %.4f*c \n",modv,modv/c);
	modB=Math->Mod(B0);
	printf("B0=%.3e T = %.3e G\n",modB,modB/1E-4);
	
	rL0=Phys->Larmor(rL,m_0,q,modv,modB);
	printf("rL=%.3e m = %.3e UA\n\n",rL0,rL0*kp);

	//Min and max for dipole and relativistc propagation

//	x_max=r0[0]+50*rL0,y_max=r0[1]+50*rL0,z_max=r0[2]+50*rL0;
//	x_min=r0[0]-50*rL0,y_min=r0[1]-50*rL0,z_min=r0[2]-0*rL0;
	
//	x_max=O[0]+10*rL0,y_max=O[1]+10*rL0,z_max=O[2]+10*rL0;
//	x_min=O[0]-10*rL0,y_min=O[1]-10*rL0,z_min=O[2]-0*rL0;

	//Min and max for the galatic field
	
	x_max=r0[0]+0.3*rs,y_max=r0[1]+0.3*rs,z_max=r0[2]+0.3*rs;
	x_min=r0[0]-0.3*rs,y_min=r0[1]-0.3*rs,z_min=r0[2]-0*rs;

	r[3]=z_min;
		
//	printf("UA=%.3e \n",UA);

    fprintf(outfile,"%.3e %.3e %.3e %.3e %.3e %.3e\n",x_min*kp,y_min*kp,z_min*kp,x_max*kp,y_max*kp,z_max*kp); //Give the max and minimum values
//   fprintf(outfile,"%.3e %.3e %.3e %.3e %.3e %.3e\n",x_min/UA,y_min/UA,z_min/UA,x_max/UA,y_max/UA,z_max/UA); //Give the max and minimum values
    
		
    w=Phys->Angvelocity(w,q,modB,m_0);
    T=Phys->period(T,w);
    t=T/10; 
    
	printf("w=%e rad/s\n",w);	
    printf("T=%e s\n",T);
    printf("t=%e s\n",t);
    
	int NP=7*int(T/t); //Modificar n° de pontos
	printf("NP=%d \n\n",NP);
		
	modr=Math->Mod(r0);
	
//	LOOP PRINCIPAL

	while(modr<=rs)
//	for(i=1;i<=10;i++)
	{
		modr=Math->Mod(r0);	
		
		counter=counter+1;
		
		printf("counter=%d, r=%.4e \n",counter,modr);
		
		if(counter<=30) printf("r0x=%.6e r0y=%.6e r0z=%.6e\n",r0[0],r0[1],r0[2]);
		
//		Phys->Dip(BD,r0);  //Makes the dipolar field
		
//		Phys->VarB(B,B0,r0,rL0); //Makes the variable field 1/r
		
//		m=Phys->Relat(m,m_0,v0); //Calculate the relativistic mass

		Phys->Galactic(BG,r0);
		
		if(counter<=30) printf("BGx=%.6e, BGy=%.6e, BGz=%.6e \n",BG[0],BG[1],BG[2]);
		
//		if(i<=30) printf("m=%.4e \n",m);
		
		modB=Math->Mod(BG);
		
		Math->Projections(vpar,vper,v0,BG);
		
		if(counter<=30) printf("modB=%.6e\n",modB);
		
		modvpar=Math->Mod(vpar);
		modvper=Math->Mod(vper);
		
		//if(counter<=30)	printf("modvper=%.6e, modvpar=%.6e \n",modvper,modvpar); 
		
		rL=Phys->Larmor(rL,m,q,modvper,modB);		
		//if(i<=30)	printf("B=%.16e, rL=%.16e \n",modB,rL);
		w=Phys->Angvelocity(w,q,modB,m);
		//vol=(w*t*i)/(2*pi);
		T=Phys->period(T,w);
			
		if(counter<=30) printf("T=%.6e, w=%.6e, rL=%.6e \n",T,w,rL);
		
		fprintf(outfile2,"%.3e %.3e %.3e\n",t*counter,rL*kp,T); //tempo, raio de larmor, voltas
		
		angle=Phys->Pitch(vper,vpar);
				
		if(counter<=30) printf("angle=%.16e \n",angle);
		
		Phys->FLorentz(F,q,vper,BG); //Chama a função, passa os valores com os nomes que eu utilizo.
		
		
		for(k=0;k<=2;k++) {
			a[k]=F[k]/m;
		}
		
	//	if(i<=3)	printf("F=(%.3e,%.3e,%.3e) a=(%.3e,%.3e,%.3e) v0=(%.3e,%.3e,%.3e)\n",F[0],F[1],F[2],a[0],a[1],a[2],v0[0],v0[1],v0[2]);
		
		Phys->Velocity(v,vper,a,t); 
		
		if(counter<=30) printf("vX=%.6e vY=%.6e vZ=%.6e\n",v[0],v[1],v[2]);
		
		modv=Math->Mod(v); 
		
	// 	if(i<=30) printf(" modv=%.6e \n",modv);
//		if(i<=3)	 printf(" modv=%.6e \n",modv);
//		if(i<=3)	 printf(" v=(%.3e,%.3e,%.3e) \n",v[0],v[1],v[2]);
		
		//printf("modvpar=%.7e \n",modvpar);
	
		for(k=0;k<=2;k++) v[k]=v[k]+vpar[k];
		
//		if(i<=30) printf("vX=%.6e vY=%.6e vZ=%.6e\n",v[0],v[1],v[2]);
		
		//if(i<=30) printf("vj=%.6e \n",v);
		
		/*double vnorm[3];
		
		Math->Norm(vnorm,v); 
	
  		double v0n=Math->Mod(v0); 
		
  		for(k=0;k<=2;k++) v[k]=vnorm[k]*v0n;*/
  		
		Phys->Trajectory(r,r0,v0,a,t);
		
//		if(i<=30) printf("rX=%.6e rY=%.6e rZ=%.6e\n",r[0],r[1],r[2]);
	
//  	fprintf(outfile,"%.3e %.3e %.3e \n", r[0]/UA,r[1]/UA,r[2]/UA);
		fprintf(outfile,"%.3e %.3e %.3e \n",r[0]*kp,r[1]*kp,r[2]*kp);
		
		/*
		if(r[0]<x_min){
			x_min=r[0];
		}
		if(r[1]<y_min){
			y_min=r[1];
		}
		if(r[2]<z_min){
			z_min=r[2];
		}
		if(r[0]>x_max){
			x_max=r[0];
		}
		if(r[1]>y_max){
			y_max=r[1];
		}
		if(r[2]>z_max){
			z_max=r[2];
		}
		*/
		
		for(k=0;k<=2;k++) {
			v0[k]=v[k];
			r0[k]=r[k];
		}
		
	//	if(modr>=rs){
	//		printf("\n\n*******Reached the sun orbit*******\n\n");
	//		printf("Numero de iteradas=%d\n",counter);
	//		break;
	//	}
				
	}//fim do loop principal
	
//	printf("...\ni=%2d\n v=(%.3e,%.3e,%.3e) r=(%.3e,%.3e,%.3e) modv=%.16e \n",i,v[0],v[1],v[2],r[0],r[1],r[2],modv);
//	rewind(outfile);
//	fprintf(outfile,"%.3e %.3e %.3e %.3e %.3e %.3e\n",x_min/UA,y_min/UA,z_min/UA,x_max/UA,y_max/UA,z_max/UA); //Give the max and minimum values
	

	fclose(outfile2);
	fclose(outfile);
    printf("\n\n*******END OF RUN*******\n\n");
 	
	system("pause");
	return 0;
}
