void msmgwb2(int nm, int ns, int NGL, double temperature0, double output_parameter1[output_point_number][1+GN1+GN2],
             double Lh, double pressure1, double temperature1, double mole_concentration11, double mole_concentration12,
             double dL, double pressure2, double temperature2, double mole_concentration21, double mole_concentration22)
{
	double KA1[2][GN1][3][20],KA2[2][GN2][3][20];
	double gi[20][3];
	double b0,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17;
	int i,j,k,m,n;	
	double cutoff_wavenumber,wavenumber_interval;
	double  Lf[output_point_number];
	FILE *fp,*fp1;
		
	
	
	gka(temperature1,nm,ns,NGL,temperature0,KA1[0],KA2[0]);	
	gka(temperature2,nm+1,ns,NGL,temperature0,KA1[1],KA2[1]);	
	
	

  for(k=0;k<output_point_number;k++)
  	Lf[k]=dL*(k+0.0);  
  for(i=0;i<output_point_number;i++)
  {
  	output_parameter1[i][0]=Lf[i];
  	for(j=0;j<GN1+GN2;j++)
    	output_parameter1[i][1+j]=0.0;
  }
  
  
  
	b1=EBwavenumber(lower_wave_number,upper_wave_number,temperature1);
	b2=EBwavenumber(lower_wave_number,upper_wave_number,temperature2);
    
	for(m=0;m<gn[ns];m++)		
		for(k=0;k<NGL;k++)
		{
  		if(ns==0)
  		{
  			b3=KA1[0][m][0][k]*mole_concentration11*pressure1;  
  			b4=KA1[0][m][1][k];	
  			b5=KA1[0][m][2][k]*mole_concentration11*pressure1; 	
  			
  			b8=KA1[1][m][2][k]*mole_concentration21*pressure2;	
  		}
  		if(ns==1)
  		{
  			b3=KA2[0][m][0][k]*mole_concentration12*pressure1;  
  			b4=KA2[0][m][1][k];	
  			b5=KA2[0][m][2][k]*mole_concentration12*pressure1; 	
  			
  			b8=KA2[1][m][2][k]*mole_concentration22*pressure2;	
  		}  		     					
			
			for(i=0;i<output_point_number;i++)
			{
  			b9=b5*Lh;
  			b10=b8*Lf[i];
  			output_parameter1[i][1+m]+=b3/(b5+1.0e-12)*(E3(b10)-E3(b9+b10))*b4;
    	}
    } 
      	
}