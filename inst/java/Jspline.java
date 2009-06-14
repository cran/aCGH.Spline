
import java.io.File;
import java.io.FilenameFilter;
import java.util.Arrays;
import java.util.Vector;
import java.util.regex.*;
import java.io.*;
import java.lang.Math.*;
/**
 *
 * Author: Tomas William Fitzgerald
 */
public class Jspline {

	public static String Staticdata[][] = null;
	public static int iterate = 0;
	public static int knots = 0;
	
	public static String Outfile = "";
	public static String OutfileName = "";

public Jspline() {

}

	
public static double[] SetValues(double[] list) {
	iterate = (int)list[0];
	knots = (int)list[1];
	return(list);
}

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


public static double[] RunSpline(double[]x, double[] y) {

       double[] result = null;
       result = new double[x.length];

       int nu=x.length;
       double e = 0;
       int colum = 2;
       int len=x.length;
	   
	   // This is the number of knots
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

       // start offset - 5 iterations
       for (int p = 0; p<Offa.length;p++) {

           for (int i=0;i<ind.length;i++) {
               ind[i] = ind[i] - Offa[p];
           }

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

         for(int i = 0; i < xPol.length; i++) {
            for(int j = i+1; j < xPol.length; j++) {
                if(xPol[j] < xPol[i]) {
                    double temp = xPol[i];
                    xPol[i] = xPol[j];
                    xPol[j] = temp;
                } if (yPol[j] < yPol[i]) {
                    double temp2 = yPol[i];
                    yPol[i] = yPol[j];
                    yPol[j] = temp2;
                }

            }

         }

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
                     // count =1;
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

     /// TO RUN JAVA SPLINE
       double[] xxx = JsplineFun(uniX,uniY,xx);

       for (int i =0;i<result.length;i++) {

			result[i] += xxx[i] /iterate;

		}

	 xxx=null;
     uniX = null;
     uniY = null;
     System.gc();
     } // end offset!!

	   xx=null;
       sX = null;
       sY = null;
       System.gc();

return(result);
}

    public static double[] JsplineFun(double[] xX, double[] yY, double[] testX) {

    int sN = xX.length;
    int sNN = testX.length;
    int i, j, k, l;
    double t, ul, dx, tmp;
    int nm1 = sN - 1;
    int n_1 = sNN - 1;

    double[] sX = new double[sN+1];
    double[] sY = new double[sN+1];
    double[] sB = new double[sN+1];
    double[] sC = new double[sN+1];
    double[] sD = new double[sN+1];
    double[] sV = new double[sNN];
    double[] retArr = new double[sNN];

    for (int p =0;p<sNN;p++) {
        retArr[p] = testX[p];
        sV[p] = testX[p];
    }

    for (int p=sN;p>0;p--) {
        sX[p] = xX[p-1];
        sY[p] = yY[p-1];
        sB[p] = 0;
        sC[p] = 0;
        sD[p] = 0;
    }
    sX[0] = 0;
    sY[0] = 0;
    sB[0] = 0;
    sC[0] = 0;
    sD[0] = 0;

    Arrays.sort(sX);
    Arrays.sort(sY);

    sD[1] = sX[2] - sX[1];
    sC[2] = (sY[2] - sY[1])/sD[1];

    for( i=2 ; i<sN ; i++) {
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
    for(i=sN-2 ; i>1 ; i--) {
	sC[i] = (sC[i]-sD[i]*sC[i+1])/sB[i];
    }

    sC[1] = sC[sN] = 0.0;
    sB[sN] = (sY[sN] - sY[sN-1])/sD[sN-1] + sD[sN-1]*(sC[sN-1]+ 2.0*sC[sN]);

    for(i=1 ; i<=nm1 ; i++) {
	sB[i] = (sY[i+1]-sY[i])/sD[i] - sD[i]*(sC[i+1]+2.0*sC[i]);
	sD[i] = (sC[i+1]-sC[i])/sD[i];
	sC[i] = 3.0*sC[i];
    }

    sC[sN] = 3.0*sC[sN];
    sD[sN] = sD[nm1];

     for (int p=0;p<sN-1;p++) {
        sX[p] = sX[p+1];
         sY[p] = sY[p+1];
         sB[p] = sB[p+1];
         sC[p] = sC[p+1];
         sD[p] = sD[p+1];
    }
    sX[sN] = 0;
    sY[sN] = 0;
    sB[sN] = 0;
    sC[sN] = 0;
    sD[sN] = 0;

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
	    }
	    while(j > i+1);
	}
   // if (sX[i] !=0) {
	dx = ul - sX[i];
	tmp = (ul < sX[0]) ? 0.0 : sD[i];
	retArr[l] = sY[i] + dx*(sB[i] + dx*(sC[i] + dx*tmp));
  //  }
    }

    return(retArr);

}

    public static double[] INTERPOL2(double[] iX, double[] oX, double[] fX) {

       double[] rX = new double[iX.length];
       double[] soX = new double[oX.length];
       double[] sfX = new double[fX.length];

       for (int i=0;i<soX.length;i++) {
           soX[i] = oX[i];
           sfX[i] = fX[i];
       }
	    
		oX = null;
        fX = null;
		System.gc();
        
		Arrays.sort(soX);
        Arrays.sort(sfX);

        for (int i=0;i<iX.length;i++) {
            
			int pin = Isearch(soX,iX[i]);
            rX[i] = iX[i] - soX[pin] + sfX[pin];
        }
		
        iX = null;
        sfX = null;
        soX = null;
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
		Staticdata = new String[size][7];
		
		if (rawOrnot.equals("T")) {
		rawOrnot = "TRUE";
		} else if (rawOrnot.equals("F")) {
		rawOrnot = "FALSE";
		}
		  
		DataInputStream dis = null; 
        String record = null; 
        int recCount = 0; 
        int red_int =0;
        int green_int =0;
		int proR =0;
		int proG =0;
        int loc_int =0;
		int gf_non =0;
		int rf_non =0;
		int gb_non =0;
		int rb_non=0;
		
        Pattern p = Pattern.compile("chr\\w\\w?:\\d+-\\d");
        int pin =0;
			
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
					
					  Staticdata[pin][2] = jchr;
					  Staticdata[pin][3] = start;
					  Staticdata[pin][4] = stop;
					  Staticdata[pin][5] = Integer.toString(recCount - 1);
					  Staticdata[pin][6] = "0";
					  
					  int check1 = Integer.parseInt(data[gf_non]);
					  int check2 = Integer.parseInt(data[rf_non]);
			          int check3 = Integer.parseInt(data[gb_non]);
			          int check4 = Integer.parseInt(data[rb_non]);
					 
					  if (check1 == 1 || check2 == 1 || check3 == 1 || check4 == 1) {
						Staticdata[pin][0] = "NA";
						Staticdata[pin][1] = "NA";
						Staticdata[pin][6] = "1";
					  } else {
						
							if (rawOrnot.equals("TRUE")) {
								Staticdata[pin][0] = data[red_int];
								Staticdata[pin][1] = data[green_int];
								Staticdata[pin][6] = "0";
							} else if (rawOrnot.equals("FALSE")) {
								Staticdata[pin][0] = data[proR];
								Staticdata[pin][1] = data[proG];
								Staticdata[pin][6] = "0";
							} 
					}
					  
					  if (rawOrnot.equals("TRUE")) {
					  
					  if (!Staticdata[pin][0].equals("NA") && !Staticdata[pin][2].equals("NA")) {
					  
							if (Double.parseDouble(Staticdata[pin][0]) + Double.parseDouble(Staticdata[pin][1]) < 100) {
							Staticdata[pin][0] = data[red_int];
							Staticdata[pin][1] = data[green_int];
							Staticdata[pin][6] = "1";
							}
					  
							if (Double.parseDouble(Staticdata[pin][0]) < 1) {
							Staticdata[pin][0] = "1";
							} 
							
							if (Double.parseDouble(Staticdata[pin][1]) < 1) {
							Staticdata[pin][1] = "1";
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

public static String[] ReturnHelper(String index) {
	String[] retArr = new String[Staticdata.length];
	
	int index1 = 0;
	
	if (index.equals("red")) 
		index1 = 0;
	else if (index.equals("green")) 
		index1 = 1;
	else if (index.equals("chr")) 
		index1 = 2;
	else if (index.equals("start")) 
		index1 = 3;
	else if (index.equals("stop")) 
		index1 = 4;
	else if (index.equals("index")) 
		index1 = 5;
	else if (index.equals("flag")) 
		index1 = 6;
	
		for (int i =0; i<retArr.length; i++) {
		retArr[i] = Staticdata[i][index1];
		}
	
return(retArr);
}


public static double[] OutAll(double[] red, double[] green, double[] flag) {
	 
	 double[] fakeReturn = new double[1];
     // Declare output stream
     PrintWriter pw = null;
     // Create new file instance
     File output = new File(OutfileName);

    if (output.exists()) {
        System.out.println("File already exists - overwriting!!!");
        //System.exit(0);
    }

     DataInputStream dis = null;
        String record = null;
        int recCount = 0;
        int probe_int = 0;
        int red_int =0;
        int green_int =0;
        int loc_int =0;
        int ratio_int =0;
		int gf_non =0;
		int rf_non =0;
		int gb_non =0;
		int rb_non =0;
		int indexP =0;
		int proR =0;
		int proG =0;
		
        Pattern p = Pattern.compile("chr\\w\\w?:\\d+-\\d");
        int pin =0;

        try {

           File f = new File(Outfile);
           FileInputStream fis = new FileInputStream(f);
           BufferedInputStream bis = new BufferedInputStream(fis);
           dis = new DataInputStream(bis);
           pw = new PrintWriter(new FileOutputStream(output), true);

           while ( (record=dis.readLine()) != null ) {
               recCount++;
        // Create data output stream for input
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
                  //Matcher matc = p.matcher(data[loc_int]);
                  Matcher matc = p.matcher(record); // also changed because of checkens!!!!

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
					  double rat = red[pin] / green[pin];
					  double LogRat10 = Math.log10(rat);
					  // double LogRat2 = Math.log(Math.pow(2, correct[pin])) / Math.log(10);
                      String repl = Double.toString(LogRat10);
                      data[ratio_int] = repl; // extract the good data
					
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
                     // System.out.println(sb);
                      pin++;

                } else {
                   pw.print(record + "\n");
                }

              }

           }

        } catch (IOException e) {
           // catch io errors from FileInputStream or readLine()
           System.out.println("Uh oh, got an IOException error!" + e.getMessage());

        } finally {
           // if the files opened okay, make sure we close it
           if (dis != null) {
	      try {
                 pw.close();
                 dis.close();
	      } catch (IOException ioe) {
	      }
           }
      }
	  return(fakeReturn);
}

public static void FreeMem() {
Staticdata = null;
System.gc();
}


}
