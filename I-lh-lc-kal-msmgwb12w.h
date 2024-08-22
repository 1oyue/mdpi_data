void msmgwb4(int nm, int ns, int NGL, double temperature0, double output_parameter1[output_point_number][1+GN1+GN2], double temperaturew,
  		       double Lh, double pressure1, double temperature1, double mole_concentration11, double mole_concentration12,
  		       double Lc, double pressure2, double temperature2, double mole_concentration21, double mole_concentration22,
  		       double dL, double pressure3, double temperature3, double mole_concentration31, double mole_concentration32)
{
	double KA1[4][GN1][3][20],KA2[4][GN2][3][20];
	double gi[20][3];
	double b0,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17;
	int i,j,k,m,n;	
	double cutoff_wavenumber,wavenumber_interval;
	double  Lf[output_point_number];
	FILE *fp,*fp1;
		

	gka(temperaturew,nm,ns,NGL,temperature0,KA1[0],KA2[0]);      	
	gka(temperature1,nm+1,ns,NGL,temperature0,KA1[1],KA2[1]);	
	gka(temperature2,nm+2,ns,NGL,temperature0,KA1[2],KA2[2]);	
	gka(temperature3,nm+3,ns,NGL,temperature0,KA1[3],KA2[3]);
	


  for(k=0;k<output_point_number;k++)
  	Lf[k]=dL*(k+0.0);  
  for(i=0;i<output_point_number;i++)
  {
  	output_parameter1[i][0]=Lf[i];
  	for(j=0;j<GN1+GN2;j++)
    	output_parameter1[i][1+j]=0.0;
  }
  
  
  
  b0=EBwavenumber(lower_wave_number,upper_wave_number,temperaturew);
	b1=EBwavenumber(lower_wave_number,upper_wave_number,temperature1);
	b2=EBwavenumber(lower_wave_number,upper_wave_number,temperature2);
    
	for(m=0;m<gn[ns];m++)		
		for(k=0;k<NGL;k++)
		{
  		if(ns==0)
  		{
  			b3=KA1[0][m][0][k];  
  			b4=KA1[0][m][1][k];	
  			b5=KA1[0][m][2][k];  	
  			
  			b6=KA1[1][m][0][k]*mole_concentration11*pressure1;	   			
  			b7=KA1[1][m][1][k];	
  			b8=KA1[1][m][2][k]*mole_concentration11*pressure1;	
  			
  			b9=KA1[2][m][0][k]*mole_concentration21*pressure2;	   			
  			b10=KA1[2][m][1][k];	
  			b11=KA1[2][m][2][k]*mole_concentration21*pressure2;	    	
  			
  			b14=KA1[3][m][2][k]*mole_concentration31*pressure3;	      					
  		}
  		if(ns==1)
  		{
  			b3=KA2[0][m][0][k];  
  			b4=KA2[0][m][1][k];	
  			b5=KA2[0][m][2][k];		
  			
  			b6=KA2[1][m][0][k]*mole_concentration12*pressure1;	   			
  			b7=KA2[1][m][1][k];	
  			b8=KA2[1][m][2][k]*mole_concentration12*pressure1;	
  			
  			b9=KA2[2][m][0][k]*mole_concentration22*pressure2;	   			
  			b10=KA2[2][m][1][k];	
  			b11=KA2[2][m][2][k]*mole_concentration22*pressure2;		    	
  				
  			b14=KA2[3][m][2][k]*mole_concentration32*pressure3;	  
  		}     						
			
			for(i=0;i<output_point_number;i++)
			{
  			b15=b8*Lh;
  			b16=b11*Lc;
  			b17=b14*Lf[i];
  			output_parameter1[i][1+m]+=(b3/(b5+1.0e-12)*b4*b0*E3(b15+b16+b17)+b6/(b8+1.0e-12)*b7*b1*(E3(b16+b17)-E3(b15+b16+b17))+b9/(b11+1.0e-12)*b10*b2*(E3(b17)-E3(b16+b17)))/b1;
    	}
    }     
     	
}