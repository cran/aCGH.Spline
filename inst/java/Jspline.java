/**
 *
 *  This program is distributed in the hope that it might be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  http://www.r-project.org/Licenses/
 *
*/


import java.util.ArrayList;
import java.io.File;
import java.io.FilenameFilter;
import java.util.Arrays;
import java.util.Vector;
import java.util.regex.*;
import java.io.*;
import java.lang.Math.*;


/**
 *
 * 	Java Class to perform:
 *
 * 	- I/O of standard aCGH microarray formats (Agilent & NimbleGen).
 * 	- Robust spline interpolation of cy5 & cy3 intensity distributions.
 * 	- Segmentation of ratio values.
 *
 * 	Author:	Tomas William Fitzgerald
 * 	Date:	15/06/2009
 * 	Email:	tf2@sanger.ac.uk
 *
 * 	Modified(tf2):	11/01/2010
 * 
 *	Note - This class is written specically to interact with the R package, aCGH.Spline.
 *
 *
**/


public class Jspline {

	/************** Spline params  *************/
	public static String Staticdata[][] = null;
	public static int iterate = 0;
	public static int knots = 0;
	/******************************************/
	
	/********** Segmentation params ***********/
 	public static double tr1 = 2;
    public static int len1 = 25;
    public static int len2 = 5;    
    public static int wlen1 = 100;
    public static int wlen2 = 500;
    public static double tr2 = 2;
    public static int wlen3 = 10;   
    public static int siZ = 500;
	/******************************************/
	
	/************* Output params **************/
	public static String Outfile = "";
	public static String OutfileName = "";
	/******************************************/

	/* Default constructor */
	public Jspline() {

	}

	/* Set spline params */
	public static double[] SetValues(double[] list) {
		iterate = (int)list[0];
		knots = (int)list[1];
	return(list);
	}

	/* Set output params */
	public static void SetOutput(String file) {
	
	Outfile = file;
	String justFileName;
	String path;
	
	int index = file.lastIndexOf('/');

        if (index >= 0) {
            justFileName = file.substring(index + 1);
            path = file.substring(0, index) + "/";
			String name = (String)justFileName.subSequence(1, justFileName.length());
			OutfileName = path + "S" + name;
        } else {
            index = file.lastIndexOf('\\');
            if (index >= 0) {
                justFileName = file.substring(index + 1) + "\\";
                path = file.substring(0, index);
				String name = (String)justFileName.subSequence(1, justFileName.length());
				OutfileName = path + "S" + name;
            } else {
                justFileName = file;
                path = file;
				String name = (String)justFileName.subSequence(1, justFileName.length());
				OutfileName = "S" + name; 
            }
        }
	
	}

	/* Run spline - define unique knots & iterations */ 
	public static double[] RunSpline(double[]x, double[] y) {

       double[] result = null;
       result = new double[x.length];
       int nu=x.length;
       double e = 0;
       int colum = 2;
       int len=x.length;
	   
	   // knots
       double number = Math.round(len/knots);
       number = Math.max(number, 100);
       double[] pin = new double[(int)number];

       for (int i=0;i<pin.length;i++) {
       		pin[i] = (i+1) / number;
       }
       int[] ind = new int[(int)number];

       for (int i=0;i<ind.length;i++) {
       		ind[i] = (int) (pin[i] * len);
       }

	   // offsets
       int Off = Math.round(ind[0]/iterate);
       int[] Offa = new int[iterate];
       Offa[0] = 0;
       
       for (int i=1;i<iterate;i++) {
       		Offa[i] = Off;
       }
	   double[] xx = new double[x.length];

	   for (int i=0 ; i<x.length;i++) {
	   		if (x[i] < 1)  {
				xx[i] = 1;
			} else {
				xx[i] = x[i];
			}
       }
       double[] sY = new double[x.length];

		 // geometric mean
         for (int u = 0; u<sY.length; u++) {       	
         	if (x[u]<1) 
                 x[u] =1;
            if (y[u]<1) 
                 y[u] =1;
             sY[u] = Math.exp((Math.log(x[u]) + Math.log(y[u])) / 2);
			 if (sY[u] < 1)
				sY[u] = 1;
         }

       double[] sX = new double[x.length];
       for (int i=0;i<sX.length;i++) {
           sX[i] = xx[i];
       }
       Arrays.sort(sX);
       Arrays.sort(sY);

       // start offset 
       for (int p = 0; p<Offa.length;p++) {
           
           for (int i=0;i<ind.length;i++) {
               ind[i] = ind[i] - Offa[p];
           }

	   // unique values
       double[] xPol = new double[(int)number];
       for (int i=0;i<xPol.length;i++) {
			xPol[i] = sX[(int)ind[i]-1];
       }
	   xPol[xPol.length-1] = sX[len-1];
       
       double[] yPol = new double[(int)number];	   
	   for (int i=0;i<yPol.length;i++) {
			yPol[i] = sY[(int)ind[i]-1];
	   }
	   yPol[yPol.length-1] = sY[len-1];
	   Arrays.sort(xPol);
       Arrays.sort(yPol);
		
         Vector<Double> vector1 = new Vector();
         Vector<Double> vector2 = new Vector();
         int count=0;
         int pin2=0;
         double tempA=0;
         double tempB=0;

           for (int i=0;i<xPol.length;i++) {
               double tempX =xPol[i];
               double tempY =yPol[i];

                if (i == xPol.length-1) {
                  double ave1 = 0;
                  double ave2=0;
                  tempX = xPol[i];
                  tempY =yPol[i];
                  if (count == 0) {
                      count =1;
                      ave1 = tempX;
                      ave2 = tempY;
                  } else if (count > 0) {
                  	ave1 = tempA / count;
                  	ave2 = tempB / count;
                  }
                 vector2.add(pin2, ave1);
                 vector1.add(pin2, ave2);
                 tempA =0;
                 count=0;

                } else if (i < xPol.length-1) {
                    int k = i+1;

                    if (xPol[i] == xPol[k]) {
                    	tempA =  tempA + tempX;
                    	tempB = tempB + tempY;
                    	count++;
                    } else {
                    	double ave1=0;
                    	double ave2=0;

                    		if (count == 0) {
                     			ave1 = tempX;
                     			ave2 = tempY;
                    		} else if(count>0) {
                    			ave1 = tempA / count;
                    			ave2 = tempB / count;
                    		}
                    		
                    	vector2.add(pin2, ave1);
                    	vector1.add(pin2, ave2);
                    	tempA =0;
                    	tempB =0;
                    	count=0;
                    	pin2++;
                    }
                }
           }

         xPol = null;
         yPol = null;
         System.gc();
         int size = vector1.size() ;
         double[] uniY = new double[size];
         double[] uniX = new double[size];

         for (int i=0;i<uniX.length;i++) {
             double temp = vector1.get(i);
             uniY[i] = temp;
             temp = vector2.get(i);
             uniX[i] = temp;
         }

		vector1 = null;
		vector2 = null;
		System.gc();

	   // Jspline
       double[] xxx = JsplineFun(uniX,uniY,xx);

	   // iterations
       for (int i =0;i<result.length;i++) {
			result[i] += xxx[i] /iterate;
		}

	 	xxx=null;
     	uniX = null;
     	uniY = null;
     	System.gc();
     	
     	} // end offset!!

	   // clean-up!! (this is needed!)
	   xx=null;
       sX = null;
       sY = null;
       System.gc();

	return(result);
	}

	/* The main spline function - see also spline.fun (R function) */
    public static double[] JsplineFun(double[] xX, double[] yY, double[] testX) {

    	int sN = xX.length;
    	int sNN = testX.length;
    	int i, j, k, l;
    	double t, ul, dx, tmp;
    	int nm1 = sN - 1;
    	int n_1 = sNN - 1;
    	
    	double[] sV= new double[sNN]; double[] retArr = new double[sNN];
    	double[] sX=new double[sN+1]; double[] sY=new double[sN+1];
    	double[] sB=new double[sN+1]; double[] sC=new double[sN+1]; double[] sD=new double[sN+1];


    	for (int p =0; p<sNN; p++) {
        	retArr[p] = testX[p];
        	sV[p] = testX[p];
    	}

    	for (int p=sN; p>0; p--) {
        	sX[p] = xX[p-1];
        	sY[p] = yY[p-1];
        	sB[p] = 0; sC[p] = 0; sD[p] = 0;
    	}
    
    	sX[0] = 0; sY[0] = 0; sB[0] = 0; sC[0] = 0; sD[0] = 0;

    	Arrays.sort(sX);
    	Arrays.sort(sY);
    	sD[1] = sX[2] - sX[1];
    	sC[2] = (sY[2] - sY[1])/sD[1];

    	for(i=2 ; i<sN ; i++) {
			sD[i] = sX[i+1] - sX[i];
			sB[i] = 2.0 * (sD[i-1] + sD[i]);
			sC[i+1] = (sY[i+1] - sY[i])/sD[i];
			sC[i] = sC[i+1] - sC[i];
    	}

    	for(i=3 ; i<sN ; i++) {
			t = sD[i-1]/sB[i-1];
			sB[i] = sB[i] - t*sD[i-1];
			sC[i] = sC[i] - t*sC[i-1];
    	}
    	sC[nm1] = sC[nm1]/sB[nm1];
    
    	for(i=sN-2; i>1; i--) {
			sC[i] = (sC[i]-sD[i]*sC[i+1])/sB[i];
    	}

    	sC[1] = sC[sN] = 0.0;
    	sB[sN] = (sY[sN] - sY[sN-1])/sD[sN-1] + sD[sN-1]*(sC[sN-1]+ 2.0*sC[sN]);

    	for(i=1; i<=nm1; i++) {
			sB[i] = (sY[i+1]-sY[i])/sD[i] - sD[i]*(sC[i+1]+2.0*sC[i]);
			sD[i] = (sC[i+1]-sC[i])/sD[i];
			sC[i] = 3.0*sC[i];
    	}

    	sC[sN] = 3.0*sC[sN];
    	sD[sN] = sD[nm1];

     	for (int p=0; p<sN-1; p++) {
        	sX[p] = sX[p+1]; sY[p] = sY[p+1];
        	sB[p] = sB[p+1]; sC[p] = sC[p+1]; sD[p] = sD[p+1];
    	}
    	sX[sN] = 0; sY[sN] = 0; sB[sN] = 0; sC[sN] = 0; sD[sN] = 0;

  		for(l = 0; l < sNN; l++) {
      		retArr[l] = sV[l];
  		}

    	i = 0;
    	for(l = 0; l < sNN; l++) {
			ul = retArr[l];
			
			if(ul < sX[i] || (i < n_1 && sX[i+1] < ul)) {
	    		i = 0;
	    		j = sN;
	    	
	    		do {
					k = (i+j)/2;
					if(ul < sX[k]) { j = k; }
					else { i = k; }
	    		} while(j > i+1);
			}

		dx = ul - sX[i];
		tmp = (ul < sX[0]) ? 0.0 : sD[i];
		retArr[l] = sY[i] + dx*(sB[i] + dx*(sC[i] + dx*tmp));
    	}

    return(retArr);
	}

    public static double[] INTERPOL2(double[] iX, double[] oX, double[] fX) {

        double[] rX = new double[iX.length];   
		Arrays.sort(oX);
        Arrays.sort(fX);

        for (int i=0;i<iX.length;i++) {
			int pin = Isearch(oX,iX[i]);
            rX[i] = iX[i] - oX[pin] + fX[pin];
        }
		
        iX = null;
        System.gc();

    return (rX);
    }

    public static int Isearch(double[] data, double point) {
        int ipin = 0;
	    int i = 0;
        int pin = 0;
		while(ipin == 0 && i < data.length-1) {
			double t = Math.abs(data[i] - point);
			double tt = Math.abs(data[i+1]-point);
				if (t < tt) {
					ipin = 1;
                	pin = i;
				}
            	if (i == data.length-2) {
                	ipin = 1;
                	pin = i-1;
            	}
			i++;
		}	
    return(pin);
    }
    
    /* Segmentation functions - see aCGH.RWalk for the full version (coming to CRAN soon) */
    
	public static double[] runSegN(double[] ddd) {

		// STEP 1: - crude interval estimation 
        double[] r1 = RunItApp(ddd);
        r1 = runMedian(r1, 1001);   
        double[] r2 = runMedianR(r1, 1001);
        r2 =getRunningLower(r2, 0.68);
        r1 = findMedians(ddd, r2);
        r1 = checkSmallCall(ddd,r1, r2);
        r2 = findMedians(ddd, r1);
        r2 = runMedian(r2, 1001);
        r2 = runMedianR(r2, 1001);
        
    return r2;
	}


 public static double[] RunItApp(double[] l) {

    double tr1 = 2;
    int len1 = 15;
    int len2 = 2;
    int wlen1 = 30;
    int wlen2 = 30;
    double tr2 = 1.5;
    int wlen3 = 15;

         double[] rr = new double[l.length];
         double[] v = findInterval(l);
         double[] vv = mergeInterval(v);
         double[] vvv = filterInterval(vv);
         double[] r = findMedians(l, vvv);

         rr = fineFilter(l, vvv);

    return rr;
    }

	public static double[] findInterval(double[] l) {
       
        double[] v = new double[l.length];
        double d = dLRs(l);
        double t = d*tr1;
        int e = v.length-len1;

        for (int x=0;x<e;x++) {
        	int s =0;
            v[x] =0;
            double m = Math.abs(l[x]);

            if (m>t) {    
            	for (int i=0;i<len1;i++) {
                	double m2 = Math.abs(l[x+i]);
                    if (m2>t) {
                    	s=s+1;
                    }
                }
            }	
            	if (s > len2) {
                	v[x] = 1;
            	}
        }

    return v;
    }

     public static double[] mergeInterval(double[] v) {

         double[] v1 = new double[v.length];
         double[] v2 = new double[v.length];
         double[] vv = new double[v.length];
         int e = v.length-wlen1;

         for (int x=0;x<e;x++) {
             int s =0;
                for (int i=0;i<wlen1;i++) {
                    if (v[i+x] == 1)
                        s=1;
                }
             v1[x] =s;
         }

         for (int x=e;x>=wlen1;x--) {
             int s =0;
                for (int i=wlen1-1;i>=0;i--) {
                    if (v[i+x] == 1)
                        s=1;
                }
             v2[x] =s;
         }

         for (int i=0;i<vv.length;i++) {
             vv[i] = v[i];
             int t = (int) (v1[i] + v2[i]);
             	if(t == 1) {
                 	vv[i] =0;
             	} else if (t == 2) {
                 	vv[i] =1;
             	}
         }

          for (int k=vv.length-wlen1;k>=0;k--) {
          	if (vv[k]==0) {
               for (int j=k;j<k+wlen1;j++) {
                    vv[j] =0;
                }
            }
        }

     return vv;
     }

     public static double[] filterInterval(double[] v) {

         double[] cc = new double[v.length];
         int e = v.length-wlen2;

         for (int x=wlen2;x<e;x++) {

             int s1 =1;
             int s2 =1;
             cc[x] = v[x];
             
                 for(int i=x-wlen2;i<x;i++) {
                     if (v[i] == 0) {
                         s1 = 0;
                     }
                 }
                 for (int k=x+wlen2;k>=x;k--) {
                     if (v[k] == 0) {
                         s2 = 0;
                     }
                 }

             if(s1 == 0 || s2 == 0) {
                 cc[x] =0;
             } 
         }

        for (int k=0;k<cc.length;k++) {

            if (cc[k]==1) {
                for (int j=k;j>=k-wlen2;j--) {
                    cc[j] =1;
                }
            }
        }

        for (int k=cc.length-1;k>=0;k--) {

            if (cc[k]==1) {
               for (int j=k;j<k+wlen2;j++) {
                    cc[j] =1;
                }
            }
        }

     return cc;
     }


     public static double[] findMedians(double[] l, double[] v) {

         boolean t = true;
         int k=0;
         int pin1 =0;
         int pin2 = 0;
         Vector meds = new Vector();
         double[] r = new double[l.length];
         v[0] = 0;

         do {        

                while(pin1 < v.length-1 && v[pin1] == 0) {
                    pin1++;
                }

                if (v[pin1]==1) {
                    double[] me = new double[pin1-pin2];
                    System.arraycopy(l,k , me, 0, me.length);
                    double md = median(me);
                    meds.add(md);
                    k = pin1;
                    pin2 = pin1;
                }

                 while(pin1 < v.length-1 && v[pin1] == 1) {
                    pin1++;
                }

                if (v[pin1]==0) {
                    double[] me = new double[pin1-pin2];
                    System.arraycopy(l,k, me, 0, me.length);
                    double md = median(me);
                    meds.add(md);
                    k=pin1;
                    pin2 = pin1;
                }

                if (pin1 >= v.length-1) {  t = false; break; }
                
         } while(t);

         t = true;
         k=0;
         pin1 =0;
         pin2 = 0;
         ArrayList[] MEDS = new ArrayList[meds.size()];
         for (int i=0;i<MEDS.length;i++) {
                MEDS[i] = new ArrayList();
            }

            do {

                if (pin1 >= v.length-1) {  t = false; break; }

                while(pin1 < v.length-1 && v[pin1] == 0) {
                    pin1++;
                }

                if (v[pin1]==1) {
                    for (int j=k;j<pin1;j++) {
                        r[j] = Double.parseDouble("" + meds.get(pin2));
                        MEDS[pin2].add(l[j]);
                    }
                    k=pin1;
                    pin2++;
                }

                while(pin1 < v.length-1 && v[pin1] == 1) {
                    pin1++;
                }   

                if (v[pin1]==0) {
                    for (int j=k;j<pin1;j++) {
                        r[j] = Double.parseDouble("" + meds.get(pin2));
                        MEDS[pin2].add(l[j]);
                    }
                    k=pin1;
                    pin2++;
                }           

            } while(t);

            r[r.length-1] = Double.parseDouble("" + meds.get(pin2-1));
            MEDS[pin2-1].add(l[l.length-1]);


            r = extendFilter(l, r, 20);

     return r;
     }

	public static double[] extendFilter(double[] l, double[] r, int wen) {

            for (int i=0;i<l.length-wen;i++) {
                boolean ch = true;

                if (r[i]<0) {
                    for(int j=i;j<i+wen;j++) {
                        if(l[j]<r[i]) {
                            ch = false;
                        }
                    }
                } else if (r[i]>0) {
                    for(int j=i;j<i+wen;j++) {
                        if(l[j]>r[i]) {
                            ch = false;
                        }
                    }
                }
                if (ch) {
                    double p = 2;
                    for(int kk=i;kk<r.length;kk++) {
                        if (r[i]!=r[kk]) {
                            p = r[kk];
                            break;
                        }
                    }
                    r[i] = p;
                }
            }

            for (int i=l.length-1;i>wen;i--) {
                boolean ch = true;

                if (r[i]<0) {
                    for(int j=i;j>i-wen;j--) {
                        if(l[j]<r[i]) {
                            ch = false;
                        }
                    }
                } else if (r[i]>0) {
                    for(int j=i;j>i-wen;j--) {
                        if(l[j]>r[i]) {
                            ch = false;
                        }
                    }
                }
                if (ch) {
                    double p = 2;
                    for(int kk=i;kk>0;kk--) {
                        if (r[i]!=r[kk]) {
                            p = r[kk];
                            break;
                        }
                    }
                    r[i] = p;
                }
            }

     return r;
     }

	public static double[] fineFilter(double[] l, double[] v) {

         boolean t = true;
         int k=0;
         int pin1 =0;
         int pin2 = 0;
         Vector meds = new Vector();
         v[0] = 0;

         do {
             
                while(pin1 < v.length-1 && v[pin1] == 0) {
                    pin1++;
                }

                if (v[pin1]==1) {
                    double[] me = new double[pin1-pin2];
                    System.arraycopy(l,k , me, 0, me.length);
                    double md = median(me);
                    meds.add(md);
                    k = pin1;
                    pin2 = pin1;
                }

                 while(pin1 < v.length-1 && v[pin1] == 1) {
                    pin1++;
                }
                

                if (v[pin1]==0) {
                    double[] me = new double[pin1-pin2];
                    System.arraycopy(l,k, me, 0, me.length);
                    double md = median(me);
                    meds.add(md);
                    k=pin1;
                    pin2 = pin1;
                }

                if (pin1 >= v.length-1) {  t = false; break; }
                
         } while(t);

         t = true;
         k=0;
         pin1 =0;
         pin2 = 0;
         ArrayList[] MEDS = new ArrayList[meds.size()];
         for (int i=0;i<MEDS.length;i++) {
                MEDS[i] = new ArrayList();
         }

            do {

                if (pin1 >= v.length-1) {  t = false; break; }

                while(pin1 < v.length-1 && v[pin1] == 0) {
                    pin1++;
                }

                if (v[pin1]==1) {
                    for (int j=k;j<pin1;j++) {
                        MEDS[pin2].add(l[j]);
                    }
                    k=pin1;
                    pin2++;
                }

                while(pin1 < v.length-1 && v[pin1] == 1) {
                    pin1++;
                }   

                if (v[pin1]==0) {
                    for (int j=k;j<pin1;j++) {
                        MEDS[pin2].add(l[j]);
                    }
                    k=pin1;
                    pin2++;
                }           

            } while(t);

         MEDS[pin2-1].add(l[l.length-1]);
         ArrayList[] NEW_MEDS = new ArrayList[MEDS.length];

         for (int i=0;i<MEDS.length;i++) {
                NEW_MEDS[i] = new ArrayList();
         }

          Vector pinner = new Vector();

         for (int i=0;i<MEDS.length;i++) {
              double[] temp = new double[MEDS[i].size()];
                    for (int j=0;j<MEDS[i].size();j++) {
                        temp[j] = Double.parseDouble("" + MEDS[i].get(j));;
                    }

              double d = dLRs(temp)*tr2;
              double rp68 = quantile(Absol(temp), 0.68);

                    for (int a=0;a<temp.length;a++) {
                        NEW_MEDS[i].add(temp[a]);

                        boolean check = true;

                        if (a<temp.length-wlen3) {

                            for (int tt=a+1;tt<a+wlen3;tt++) {

                                double dif = temp[a] - temp[tt];

                                if (dif<d) {
                                    check = false;
                                }
                            }
                        } else if (a>temp.length-wlen3) {

                            for (int tt=temp.length-1;tt>temp.length-wlen3;tt--) {

                                double dif = temp[a] - temp[tt];

                                if (dif<d) {
                                    check = false;
                                }
                            }

                        }
                 
                        if(check) {
                            pinner.add(NEW_MEDS[i].get(a));
                        } else {
                             pinner.add(0);
                        }

                    }
                }

         double[] medis = new double[pinner.size()];
         for (int i=0;i<pinner.size();i++) {
            medis[i] = Double.parseDouble("" + pinner.get(i));
         }

         medis[0] = 0;
         t = true;
         k=0;
         pin1 =0;
         pin2 = 0;
         double[] r = new double[l.length];
         Vector fin = new Vector();

         do {

                while(pin1 < medis.length-1 && medis[pin1] == 0) {
                    pin1++;
                }

                if (medis[pin1]!=0) {
                    double[] me = new double[pin1-pin2];
                    System.arraycopy(l,k , me, 0, me.length);
                    double md = median(me);
                    fin.add(md);
                    k = pin1;
                    pin2 = pin1;
                }

                 while(pin1 < medis.length-1 && medis[pin1] != 0) {
                    pin1++;
                }

                if (medis[pin1]==0) {
                    double[] me = new double[pin1-pin2];
                    System.arraycopy(l,k, me, 0, me.length);
                    double md = median(me);
                    fin.add(md);
                    k=pin1;
                    pin2 = pin1;
                }

                if (pin1 >= medis.length-1) {  t = false; break; }

         } while(t);

         t = true;
         k=0;
         pin1 =0;
         pin2 = 0;

            do {

                if (pin1 >= medis.length-1) {  t = false; break; }

                while(pin1 < medis.length-1 && medis[pin1] == 0) {
                    pin1++;
                }

                if (medis[pin1]!=0) {
                    for (int j=k;j<pin1;j++) {
                        r[j] = Double.parseDouble("" + fin.get(pin2));
                    }
                    k=pin1;
                    pin2++;
                }

                while(pin1 < medis.length-1 && medis[pin1] !=0) {
                    pin1++;
                }

                if (medis[pin1]==0) {
                    for (int j=k;j<pin1;j++) {
                        r[j] = Double.parseDouble("" + fin.get(pin2));
                    }
                    k=pin1;
                    pin2++;
                }

            } while(t);

            r[r.length-1] = Double.parseDouble("" + fin.get(pin2-1));
            double[] a = Absol(r);
            double d = quantile(a, 0.68);
            a = Absol(l);
            double dl = quantile(a, 0.68);

            for (int i=0;i<r.length-1;i++) {
                double diff = Math.abs(r[i] - r[i+1]);
                if (diff > d) {
                    r[i] =0;
                }                           
            }

           for(int i=1;i<r.length;i++) {
                if (r[i]==0) {
                    r[i] = r[i-1];
                }
            }

     return r;
     }

	public static double[] getRunningLower(double[] r, double per) {

         double[] d = new double[r.length];
         for(int i=0;i<r.length-1;i++) {
             d[i] = Math.abs(r[i] - r[i+1]);
         }
         d[d.length-1] = 0;

         Vector v = new Vector();
         for (int i=0;i<d.length;i++) {
             if(d[i] > 0) {
                 v.add(d[i]);
             }
         }

         double[] t = new double[v.size()];
         for (int i=0;i<t.length;i++) {
             t[i] = Double.parseDouble("" + v.get(i));
         }

         double q = quantile(t, per);
         for (int i=0;i<d.length; i++) {
             if(d[i]<q) {
                 d[i] = 0;
             } else {
                 d[i] = 1;
             }
         }

     return d;
     }

     public static double[] checkSmallCall(double[] l, double[] r, double[] rr) {

         double[] a = Absol(l);
         double t = quantile(a, 0.68)*tr1;

         for (int i=0;i<l.length;i++) {

             double d = Math.abs(l[i] - r[i]);
             if (d>t) {
                 rr[i] = 1;
             }
         }

         double[] p = new double[rr.length];
         for(int i=1;i<rr.length-1;i++) {
             p[i] = rr[i];
            if(rr[i]==1) {
                if(rr[i]==rr[i-1]) {
                    p[i] = 0;
                } if (rr[i+1]==0) {
                    p[i] = 1;
                }
            } 

         }

     return p;
     }

	
	/************** Stats functions ****************/

 	public static double[] Absol(double[] values) {      
        double[] a = new double[values.length];
        for (int i =0; i<values.length ; i++) {
            a[i] = Math.abs(values[i]);
        }  
    return a;    
    }
    
    public static float Round(float Rval, int Rpl) {
        float p = (float)Math.pow(10,Rpl);
        Rval = Rval * p;
        float tmp = Math.round(Rval);
    return (float)tmp/p;
    }

    public static double mean(double[] m) {
      double men = 0;
      int l = m.length;  
      for (int i = 0 ; i<m.length ; i++) {
      men = men + m[i];     
      }     
    return men / l;
    }

  	public static double[] runMean(double[] m, int wind) {
      double[] retM = new double[m.length];
      for (int i=0;i<m.length-wind;i++) {
          double[] temp = new double[wind];
          System.arraycopy(m, i, temp, 0, temp.length);
          double men = mean(temp);
          for (int j=i;j<i+wind;j++) {
              retM[j] = men;
          }
      }
  	return(retM);
  	}

  	public static double median(double[] values) {
      double[] copy = new double[values.length];
      System.arraycopy(values, 0, copy, 0, copy.length);
      java.util.Arrays.sort(copy); 
      int middle = copy.length/2;
    	if (copy.length%2 == 1) {
      		return copy[middle];
      	} else {
       		return (copy[middle-1] + copy[middle]) / 2.0;
    	}	
	}

  	public static double[] runMedian(double[] m, int wind) {
      double[] retM = new double[m.length];
      for (int i=0;i<m.length-wind;i++) {
          double[] temp = new double[wind];
          System.arraycopy(m, i, temp, 0, temp.length);
          double men = median(temp);
          for (int j=i;j<i+wind;j++) {
              retM[j] = men;
          }
      }
  	return(retM);
  	}

  	public static double[] runMedianR(double[] m, int wind) {
      double[] retM = new double[m.length];
      double[] rev = new double[m.length];
      int pin =0;
      for (int i=m.length-1;i>=0;i--) {
          rev[pin] = m[i];
          pin++;
      }
      for (int i=0;i<rev.length-wind;i++) {
          double[] temp = new double[wind];
          System.arraycopy(rev, i, temp, 0, temp.length);
          double men = median(temp);
          for (int j=i;j<i+wind;j++) {
              retM[j] = men;
          }
      }
      pin =0;
      for (int i=m.length-1;i>=0;i--) {
          rev[pin] = retM[i];
          pin++;
      }
  	return(rev);
  	}

   	public static double quantile(double[] values, double quantile) {
        double[] copy = new double[values.length];
        System.arraycopy(values, 0, copy, 0, copy.length);
        java.util.Arrays.sort(copy);
        int index = (int) (copy.length * quantile);
    return copy[index];
    }

   	public static double IQR(double[] values) {
       double q25 = quantile(values, 0.25);
       double q75 = quantile(values, 0.75);
    return q75 - q25;
   	}

   	public static double dLRs(double[] values) {
       double[] diff = new double[values.length -1];
       for (int i = 0; i<values.length -1;i++) {
           int k = i+1;
           diff[i] = values[i] - values[k];
       }
       double dLRs = IQR(diff) / 1.907745; 
   	return dLRs;
   	}


	/********************* I/O Functions ********************/

	public static int ArraySizer(String filename) {
     Pattern p = Pattern.compile("chr\\w\\w?:\\d+-\\d");
     filename.trim();
     String inLine;
     String outline = "";
     BufferedReader infile = null;
     int len=0;
        
     try {
         infile = new BufferedReader(new FileReader (filename));
         
		 while ((inLine=infile.readLine()) != null) {
             Matcher m = p.matcher(inLine);
             
			 if (m.find()) 
             len++;
         }
     }
     catch (FileNotFoundException ex) {
         System.out.println("File not found: " + filename);
     }
     catch (IOException ex) {
         System.out.println(ex.getMessage());
     }
     finally {
         try {
             if (infile != null) infile.close();
         }
         catch (IOException ex) {
             System.out.println(ex.getMessage());
         }
     }

    return len;      
    }


	public static int InAll(String file, String rawOrnot) {
		
		int size = ArraySizer(file);
		Staticdata = new String[size][8];
		
		if (rawOrnot.equals("T")) {
		rawOrnot = "TRUE";
		} else if (rawOrnot.equals("F")) {
		rawOrnot = "FALSE";
		}
		  
		DataInputStream dis = null; 
        String record = null; 
        int recCount=0; int pin = 0; 
        int red_int=0; int green_int=0; 
        int proR=0; int proG=0; int loc_int=0;
		int gf_non=0; int rf_non=0; int gb_non=0; int rb_non=0;	
        Pattern p = Pattern.compile("chr\\w\\w?:\\d+-\\d");
			
        try { 

           File f = new File(file); 
           FileInputStream fis = new FileInputStream(f); 
           BufferedInputStream bis = new BufferedInputStream(fis); 
           dis = new DataInputStream(bis);  

           while ( (record=dis.readLine()) != null ) { 
              recCount++; 
              if (recCount == 10){
                  String[] seq = record.split("\t");
              for (int i =0; i<seq.length; i++ ) {
              String temp = seq[i];
                 if (temp.equals("rBGSubSignal"))
                    red_int = i;
                 else if (temp.equals("gBGSubSignal"))
                    green_int = i;
				 else if (temp.equals("rProcessedSignal"))
					proR = i;	
				 else if (temp.equals("gProcessedSignal"))
					proG = i;				
                 else if (temp.equals("SystematicName"))
                    loc_int = i;
				 else if (temp.equals("gIsFeatNonUnifOL"))
                    gf_non = i;
				 else if (temp.equals("rIsFeatNonUnifOL"))
                    rf_non = i;
				 else if (temp.equals("gIsBGNonUnifOL"))
                    gb_non = i;
				 else if (temp.equals("rIsBGNonUnifOL"))
                    rb_non = i;             
              }
           }
              if (recCount > 10) {
                  String[] data = record.split("\t");
				  Matcher matc= p.matcher(record); // Changed because of the chickens..
                  
                  if (matc.find()) {
				  
					  String SyName = data[loc_int];

                      String[] SyArray = SyName.split(":");
					  String chr = SyArray[0];
					  String jchr = chr.substring(3,chr.length());
					  
					  if (jchr.equals("X")) {
					  jchr = "23";
					  } else if (jchr.equals("Y")) {
					  jchr = "24";
					  }
					  
                      String SyTemp = SyArray[1];
                      String[] SyLocat = SyTemp.split("-");
					  String start = SyLocat[0];
					  String stop = SyLocat[1];
					
					  Staticdata[pin][3] = jchr;
					  Staticdata[pin][4] = start;
					  Staticdata[pin][5] = stop;
					  Staticdata[pin][6] = Integer.toString(recCount - 1);
					  Staticdata[pin][7] = "0";
					  
					  int check1 = Integer.parseInt(data[gf_non]);
					  int check2 = Integer.parseInt(data[rf_non]);
			          int check3 = Integer.parseInt(data[gb_non]);
			          int check4 = Integer.parseInt(data[rb_non]);
					 
					  if (check1 == 1 || check2 == 1 || check3 == 1 || check4 == 1) {
					    Staticdata[pin][0] = "NA";
						Staticdata[pin][1] = "NA";
						Staticdata[pin][2] = "NA";
						Staticdata[pin][7] = "1";
					  } else {
						
							if (rawOrnot.equals("TRUE")) {
								double red = Double.parseDouble(data[red_int]);
					  			double green= Double.parseDouble(data[green_int]);
								double LogRat2 = Math.log(red/ green) / Math.log(2);
							    Staticdata[pin][0] = "" + LogRat2;
								Staticdata[pin][1] = data[red_int];
								Staticdata[pin][2] = data[green_int];
								Staticdata[pin][7] = "0";
							} else if (rawOrnot.equals("FALSE")) {
							    double red = Double.parseDouble(data[proR]);
					  			double green= Double.parseDouble(data[proG]);
								double LogRat2 = Math.log(red/ green) / Math.log(2);
							    Staticdata[pin][0] = "" + LogRat2;
								Staticdata[pin][1] = data[proR];
								Staticdata[pin][2] = data[proG];
								Staticdata[pin][7] = "0";
							} 
					}
					  
					  if (rawOrnot.equals("TRUE")) {
					  
					  if (!Staticdata[pin][1].equals("NA") && !Staticdata[pin][2].equals("NA")) {
					  
							if (Double.parseDouble(Staticdata[pin][1]) + Double.parseDouble(Staticdata[pin][2]) < 100) {
								double red = Double.parseDouble(data[red_int]);
					  			double green= Double.parseDouble(data[green_int]);
								double LogRat2 = Math.log(red/ green) / Math.log(2);
								Staticdata[pin][0] = "" + LogRat2;
								Staticdata[pin][1] = data[red_int];
								Staticdata[pin][2] = data[green_int];
								Staticdata[pin][7] = "1";
							}
					  
							if (Double.parseDouble(Staticdata[pin][1]) < 1) {
								Staticdata[pin][0] = "0";
								Staticdata[pin][1] = "1";
								Staticdata[pin][2] = "1";
								Staticdata[pin][7] = "1";
							} 
							
							if (Double.parseDouble(Staticdata[pin][2]) < 1) {
								Staticdata[pin][0] = "0";
								Staticdata[pin][1] = "1";
								Staticdata[pin][2] = "1";
								Staticdata[pin][7] = "1";
							}
					  }
					} 
                  pin++;
                }
                  
              }
            
           }

           
        } catch (IOException e) {        
           System.out.println("Uh oh, got an IOException error!" + e.getMessage()); 
        } finally { 
         
		if (dis != null) { 
	      try {
                 dis.close(); 
	      } catch (IOException ioe) {
	      }
		} 
        }
	return(size);
    }	  

	public static double[] OutAll(double[] log, double[] red, double[] green, double[] flag) {
	 
	 double[] fR = new double[1];
     PrintWriter pw = null;
     File output = new File(OutfileName);

    	if (output.exists()) {
        	System.out.println("File already exists - overwriting!!!");
    	}

        DataInputStream dis = null;
        String record = null;
        int recCount=0; int pin = 0;
        int probe_int=0; int red_int=0; int green_int=0; int loc_int=0; int ratio_int=0;
        int gf_non=0; int rf_non=0; int gb_non=0; int rb_non =0;
        int indexP=0; int proR=0; int proG=0;	
        Pattern p = Pattern.compile("chr\\w\\w?:\\d+-\\d");

        try {

           File f = new File(Outfile);
           FileInputStream fis = new FileInputStream(f);
           BufferedInputStream bis = new BufferedInputStream(fis);
           dis = new DataInputStream(bis);
           pw = new PrintWriter(new FileOutputStream(output), true);

        while ( (record=dis.readLine()) != null ) {
        	recCount++;
        	
        if (recCount < 10) {
        pw.print(record + "\n");
        }

        else if (recCount == 10){
            
			pw.print(record + "\n");
			String[] seq = record.split("\t");
              
			  for (int i =0; i<seq.length; i++ ) {
              String temp = seq[i];
                 
				 if (temp.equals("ProbeName"))
                    probe_int = i;
                 else if (temp.equals("LogRatio"))
                    ratio_int = i;
				 else if (temp.equals("PValueLogRatio"))
					indexP = i;					
				 else if (temp.equals("rProcessedSignal"))
					proR = i;	
				 else if (temp.equals("gProcessedSignal"))
					proG = i;									
                 else if (temp.equals("rBGSubSignal"))
                    red_int = i;
                 else if (temp.equals("gBGSubSignal"))
                    green_int = i;
                 else if (temp.equals("SystematicName"))
                    loc_int = i;
				 else if (temp.equals("gIsFeatNonUnifOL"))
                    gf_non = i;
				 else if (temp.equals("rIsFeatNonUnifOL"))
                    rf_non = i;
				 else if (temp.equals("gIsBGNonUnifOL"))
                    gb_non = i;
				 else if (temp.equals("rIsBGNonUnifOL"))
                    rb_non = i;

              }
           }

              else if (recCount > 10) {
                  String[] data = record.split("\t");
                  Matcher matc = p.matcher(record); 

                  if (matc.find()) {
                    
					if (flag[pin] == 1) {
					data[indexP]="1";
					data[gf_non]="1";
					data[rf_non]="1";
					data[gb_non]="1";
					data[rb_non]="1";
					} 
					
					if (red[pin]<1) {
					red[pin]=1;
					}
					
					if (green[pin]<1) {
					green[pin]=1;
					}
					  					                    
					  data[proR] = Double.toString(red[pin]);
					  data[proG] = Double.toString(green[pin]);
					  String rat = Double.toString(log[pin]);;
					  double log10 = 0;
					  
					  if (rat.equals("NA")) {
					  	log10 = 0;
					  	data[indexP]="1";
						data[gf_non]="1";
						data[rf_non]="1";
						data[gb_non]="1";
						data[rb_non]="1";
					  } else {
					  		double log2 = log[pin];
					  		log10 = Math.log(Math.pow(2, log2)) / Math.log(10);
					  }
								  
                      String repl = Double.toString(log10);
                      data[ratio_int] = repl;					
                      StringBuilder sb = new StringBuilder(512); // should be enough...
                      boolean first = true;

                      for(CharSequence s : data) {
                        if(first) {
                            first = false;
                        } else {
                          sb.append("\t");
                        }
                        sb.append(s);
                      }

                      pw.print(sb + "\n");
                      pin++;

                } else {
                   pw.print(record + "\n");
                }

              }

           }

        } catch (IOException e) {
           System.out.println("Uh oh, got an IOException error!" + e.getMessage());

        } finally {
           if (dis != null) {
	      try {
                 pw.close();
                 dis.close();
	      } catch (IOException ioe) {
	      }
           }
      }
	  return(fR);
	}
	
	
	
	public static int InAllNimblegen(String file, String rawOrnot, String index) {

		int size = ArraySizer(file);
		Staticdata = new String[size][8];
		
		if (rawOrnot.equals("T")) {
		rawOrnot = "TRUE";
		} else if (rawOrnot.equals("F")) {
		rawOrnot = "FALSE";
		}

            DataInputStream dis = null;
            String record = null;
            int recCount=0; int pin=0;
            int chr_int=0; int pos_int = 0;  
            int exp_int=0; int ref_int=0; int Nexp_int=0; int Nref_int=0;          
            Pattern p = Pattern.compile("chr\\w\\w?:\\d+-\\d");
            String da = "";
            double[][] ret_vec = new double[size][2];

        try {

           File f = new File(file);
           FileInputStream fis = new FileInputStream(f);
           BufferedInputStream bis = new BufferedInputStream(fis);
           dis = new DataInputStream(bis);

           while ( (record=dis.readLine()) != null ) {
              recCount++;
              if (recCount == 1){
                  String[] seq = record.split("\t");
              for (int i =0; i<seq.length; i++ ) {
                  String temp = seq[i];
                  if (temp.equals("EXP_SPATIAL"))
                    exp_int = i;
                  else if (temp.equals("REF_SPATIAL"))
                    ref_int = i;
                  else if (temp.equals("EXP_NORM"))
                  	Nexp_int = i;
                  else if (temp.equals("REF_NORM")) 
                  	Nref_int = i;
                  else if (temp.equals("CHROMOSOME")) 
					chr_int = i;
				  else if (temp.equals("CHR_POSITION")) 
					pos_int = i;		  
              }
           }
              if (recCount > 1) {
                  String[] data = record.split("\t");
                  Matcher matc = p.matcher(record);

                  if (matc.find()) {
                  
					  String chr = data[chr_int];
					  String jchr = chr.substring(3,chr.length());

					  if (jchr.equals("X")) {
					  jchr = "23";
					  } else if (jchr.equals("Y")) {
					  jchr = "24";
					  }
					  
					  String pos = data[pos_int];
					  Staticdata[pin][3] = jchr;
					  Staticdata[pin][4] = pos;
					  Staticdata[pin][5] = pos;
					  Staticdata[pin][6] = Integer.toString(recCount - 1);					
					  Staticdata[pin][7] = "0";
					 
					     if (rawOrnot.equals("TRUE")) {
								double red = Double.parseDouble(data[exp_int]);
					  			double green= Double.parseDouble(data[ref_int]);
					  			red = Math.pow(2, red);
                      			green = Math.pow(2, green);
								double LogRat2 = Math.log(red/ green) / Math.log(2);
							    Staticdata[pin][0] = "" + LogRat2;
								Staticdata[pin][1] = "" + red;
								Staticdata[pin][2] = "" + green;
								Staticdata[pin][7] = "0";
							} else if (rawOrnot.equals("FALSE")) {
							    double red = Double.parseDouble(data[Nexp_int]);
					  			double green= Double.parseDouble(data[Nref_int]);
					  			red = Math.pow(2, red);
                      			green = Math.pow(2, green);
								double LogRat2 = Math.log(red/ green) / Math.log(2);
							    Staticdata[pin][0] = "" + LogRat2;
								Staticdata[pin][1] = "" + red;
								Staticdata[pin][2] = "" + green;
								Staticdata[pin][7] = "0";
							} 
					  
					  if (rawOrnot.equals("TRUE")) {
					  
					  if (!Staticdata[pin][1].equals("NA") && !Staticdata[pin][2].equals("NA")) {
					  
							if (Double.parseDouble(Staticdata[pin][1]) + Double.parseDouble(Staticdata[pin][2]) < 100) {
								double red = Double.parseDouble(data[exp_int]);
					  			double green= Double.parseDouble(data[ref_int]);
					  			red = Math.pow(2, red);
                      			green = Math.pow(2, green);
								double LogRat2 = Math.log(red/ green) / Math.log(2);
								Staticdata[pin][0] = "" + LogRat2;
								Staticdata[pin][1] = "" + red;
								Staticdata[pin][2] = "" + green;
								Staticdata[pin][7] = "1";
							}
					  }
					  
							if (Double.parseDouble(Staticdata[pin][1]) < 1) {
								Staticdata[pin][0] = "0";
								Staticdata[pin][1] = "1";
								Staticdata[pin][2] = "1";
								Staticdata[pin][7] = "1";
							} 
							
							if (Double.parseDouble(Staticdata[pin][2]) < 1) {
								Staticdata[pin][0] = "0";
								Staticdata[pin][1] = "1";
								Staticdata[pin][2] = "1";
								Staticdata[pin][7] = "1";
							}
						}
                  pin++;
                }
                  
              }
            
           }                             

        } catch (IOException e) {
           System.out.println("Uh oh, got an IOException error!" + e.getMessage());

        } finally {
           if (dis != null) {
	      try {
                 dis.close();
	      } catch (IOException ioe) {
	      }
           }
        }
        return size;
      }

	public static double[] OutAllNimbelgen(double[] log, double[] red, double[] green, double[] flag) {
	 
	 double[] fR = new double[1];
     PrintWriter pw = null;
     File output = new File(OutfileName);

    	if (output.exists()) {
        	System.out.println("File already exists - overwriting!!!");
    	}

     	DataInputStream dis = null;
        String record = null;
        int recCount=0; int pin=0;
        int Nexp_int=0; int Nref_int=0; int chr_int=0; int pos_int=0; int ratio_int = 0;
		
        Pattern p = Pattern.compile("chr\\w\\w?:\\d+-\\d");

        try {

           File f = new File(Outfile);
           FileInputStream fis = new FileInputStream(f);
           BufferedInputStream bis = new BufferedInputStream(fis);
           dis = new DataInputStream(bis);
           pw = new PrintWriter(new FileOutputStream(output), true);

           while ( (record=dis.readLine()) != null ) {
               recCount++;

        	if (recCount < 10) {
        		pw.print(record + "\n");
        	}

        	else if (recCount == 1){
            
			pw.print(record + "\n");
			String[] seq = record.split("\t");
              
			  for (int i =0; i<seq.length; i++ ) {
              String temp = seq[i];

                  if (temp.equals("EXP_NORM"))
                  	Nexp_int = i;
                  else if (temp.equals("REF_NORM")) 
                  	Nref_int = i;
                  else if (temp.equals("CHROMOSOME")) 
					pos_int = i;
				  else if (temp.equals("POSITION")) 
					pos_int = i;					
				  else if (temp.equals("RATIO_CORRECTED"))
				  	ratio_int = i;

              }
           }

              else if (recCount > 1) {
                  String[] data = record.split("\t");
                  Matcher matc = p.matcher(record); 

                  if (matc.find()) {
                    
                    double logOut = 0;                   					
					
					if (red[pin]<1) {
						red[pin]=1;
					}
					
					if (green[pin]<1) {
						green[pin]=1;
					}
					  					                    
					  data[Nexp_int] = Double.toString(red[pin]);
					  data[Nref_int] = Double.toString(green[pin]);
					  String rat = Double.toString(log[pin]);;
					  
					  if (rat.equals("NA")) {
					  	logOut = 0;
					  } else if (flag[pin] == 0) {
					  	logOut = 0;
					  } else {
					  		logOut = log[pin];
					  }
								  
                      String repl = Double.toString(logOut);
                      data[ratio_int] = repl;					
                      StringBuilder sb = new StringBuilder(512); // should be enough...
                      boolean first = true;

                      for(CharSequence s : data) {
                        if(first) {
                            first = false;
                        } else {
                          sb.append("\t");
                        }
                        sb.append(s);
                      }

                      pw.print(sb + "\n");
                      pin++;

                } else {
                   pw.print(record + "\n");
                }

              }

           }

        } catch (IOException e) {
           System.out.println("Uh oh, got an IOException error!" + e.getMessage());

        } finally {
           if (dis != null) {
	      try {
                 pw.close();
                 dis.close();
	      } catch (IOException ioe) {
	      }
           }
      }
	  return(fR);
	}
	
	public static String[] ReturnHelper(String index) {
	String[] retArr = new String[Staticdata.length];
	
	int index1 = 0;
	
	if (index.equals("log")) 
		index1 = 0;
	else if (index.equals("red")) 
		index1 = 1;
	else if (index.equals("green")) 
		index1 = 2;
	else if (index.equals("chr")) 
		index1 = 3;
	else if (index.equals("start")) 
		index1 = 4;
	else if (index.equals("stop")) 
		index1 = 5;
	else if (index.equals("index")) 
		index1 = 6;
	else if (index.equals("flag")) 
		index1 = 7;
	
		for (int i =0; i<retArr.length; i++) {
		retArr[i] = Staticdata[i][index1];
		}
	
	return(retArr);
	}


	public static void FreeMem() {
		Staticdata = null;
		System.gc();
	}

}
