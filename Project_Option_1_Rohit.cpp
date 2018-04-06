//
//  MRC.cpp
//  ImageRead
//
//  Created by salim sazzad on 12/15/16.
//  Copyright (c) 2016 salim sazzad. All rights reserved.
//


#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <vector>
#include <cmath>
#include "MRC.h"


#define DENSITY_INFINITY 1000000
int main()
{
  MRC MRCobject;
  string InputFileName;
  cout << "Enter the name of input file (Ex. 5030_BR.mrc) " << endl;
  cin >> InputFileName;
  MRCobject.readMRCFile(InputFileName);
  MRCobject.writeDataInMRCFormat("output_File_Rohit.mrc");
  MRCobject.minDensity();
  MRCobject.maxDensity();
  MRCobject.meanDensity();
  MRCobject.showMRCFileDescription();
  return 0;
}

/*
 * Read from MRC file
 * Provide MRC file path name as an argument
 * MRC file is a binary file, so you need to read byte by byte
 * MRC file format: http://bio3d.colorado.edu/imod/doc/mrc_format.txt
*/

void MRC::readMRCFile(string inputMRCFile)
{
    printf("\nReading MRC file: %s\n", inputMRCFile.c_str());

    fstream inputFile(inputMRCFile.c_str(), ios::in|ios::binary);
    if(!inputFile)
    {
        cout<<"Can not open binary input file "<<inputMRCFile<<endl;
        exit(0);
    }

    inputFile.read(reinterpret_cast<char *>(&nx), sizeof(nx));  //number of columns of the complete image
    inputFile.read(reinterpret_cast<char *>(&ny), sizeof(ny));  //number of rows
    inputFile.read(reinterpret_cast<char *>(&nz), sizeof(nz));  //number of sections
    inputFile.read(reinterpret_cast<char *>(&mode), sizeof(mode));  //types of pixels
    inputFile.read(reinterpret_cast<char *>(&nxStart), sizeof(nxStart));  // starting point of subimage
    inputFile.read(reinterpret_cast<char *>(&nyStart), sizeof(nyStart));
    inputFile.read(reinterpret_cast<char *>(&nzStart), sizeof(nzStart));
    inputFile.read(reinterpret_cast<char *>(&mx), sizeof(mx));  // grid size in x,y,z
    inputFile.read(reinterpret_cast<char *>(&my), sizeof(my));
    inputFile.read(reinterpret_cast<char *>(&mz), sizeof(mz));
    inputFile.read(reinterpret_cast<char *>(&xLength), sizeof(xLength));
    inputFile.read(reinterpret_cast<char *>(&yLength), sizeof(yLength));  // cell sizing
    inputFile.read(reinterpret_cast<char *>(&zLength), sizeof(zLength));
    inputFile.read(reinterpret_cast<char *>(&alpha), sizeof(alpha));
    inputFile.read(reinterpret_cast<char *>(&beta), sizeof(beta));  //cell angles
    inputFile.read(reinterpret_cast<char *>(&gamma), sizeof(gamma));
    inputFile.read(reinterpret_cast<char *>(&mapc), sizeof(mapc));
    inputFile.read(reinterpret_cast<char *>(&mapr), sizeof(mapr));  //mapping
    inputFile.read(reinterpret_cast<char *>(&maps), sizeof(maps));
    inputFile.read(reinterpret_cast<char *>(&dMin), sizeof(dMin));
    inputFile.read(reinterpret_cast<char *>(&dMax), sizeof(dMax));  //scaling of data
    inputFile.read(reinterpret_cast<char *>(&dMean), sizeof(dMean));
    inputFile.read(reinterpret_cast<char *>(&ispg), sizeof(ispg));  //space group number
    inputFile.read(reinterpret_cast<char *>(&nsymbt), sizeof(nsymbt)); //bytes used for symmentry data
    inputFile.read(reinterpret_cast<char *>(extra), 25*sizeof(int));  //extra space
    inputFile.read(reinterpret_cast<char *>(&xOrigin), sizeof(xOrigin));
    inputFile.read(reinterpret_cast<char *>(&yOrigin), sizeof(yOrigin));
    inputFile.read(reinterpret_cast<char *>(&zOrigin), sizeof(zOrigin));
    inputFile.read(reinterpret_cast<char *>(map), 4*sizeof(char));
    inputFile.read(reinterpret_cast<char *>(&machineStamp), sizeof(machineStamp));
    inputFile.read(reinterpret_cast<char *>(&rms), sizeof(rms));
    inputFile.read(reinterpret_cast<char *>(&nlabl), sizeof(nlabl));
    inputFile.read(reinterpret_cast<char *>(label), 10*80*sizeof(char));

    if (nx < 0 || ny < 0 || nz < 0) {
        printf("negetive nx , ny, nz value is MRC map");
        return;
    }

    cube = new float **[nx];

    for(int i = 0; i < ny; i++)
    {
        cube[i] = new float *[ny];
        for(int j = 0; j < nz; j++)
            cube[i][j] = new float [nz];
    }

    for(int i = 0; i < nz; i++){
        for(int j = 0; j < ny; j++){
            for(int k = 0; k < nx; k++)
            {
                inputFile.read(reinterpret_cast<char *>(&cube[k][j][i]), sizeof(float));

            }
        }
    }
    inputFile.close();
}


/*
 * Write data in MRC format.
 * file is the output file path.
 * set value of all the header field(nx,ny.......,label) before writting to outputMRCFile in binary Format.
 * e.g nx(number of column),ny (num of row),nz(number of section).
 * to calculate dMin, dMax, dMean you may use provided function.
*/

void MRC::writeDataInMRCFormat(string outputMRCFile)
{
    int c1, c2, c3;
    c1 = c2 = c3 = 0;
    int side, angle;

    cout << "Enter the coordinates for the center of chopping box " << endl;
    cout << "Xc: ";
    cin >> c1;
    cout << "Yc: ";
    cin >> c2;
    cout << "Zc: ";
    cin >> c3;
    cout << endl;
    cout << "Enter side S of the chopping box  ";
    cin >> side;
    cout << "Enter angle in degrees to rotate aroung z-axis  ";
    cin >> angle;
    angle = 360 - angle;

    if(c1>nx || c2>ny || c3>nz)
    {
      cout << "**** Entered coordinates are out of range ****"<< endl;
      cout << "*********** Terminating Program ***********" << endl << endl;
      return;
    }
    int nxS = c1 - side/2;
    int nxE = c1 + side/2;
    int nyS = c2 - side/2;
    int nyE = c2 + side/2;
    int nzS = c3 - side/2;
    int nzE = c3 + side/2;

    printf("\nWritting MRC data to binary file: %s\n", outputMRCFile.c_str());

    fstream outputFile(outputMRCFile.c_str(), ios::out|ios::binary);
    if(!outputFile)
    {
        cout<<"Can not open binary output file "<<outputMRCFile<<endl;
        exit(0);
    }
    nx = side;
    ny = side;
    nz = side;

    outputFile.write(reinterpret_cast<char *>(&nx), sizeof(nx));
    outputFile.write(reinterpret_cast<char *>(&ny), sizeof(ny));
    outputFile.write(reinterpret_cast<char *>(&nz), sizeof(nz));
    outputFile.write(reinterpret_cast<char *>(&mode), sizeof(mode));
    outputFile.write(reinterpret_cast<char *>(&nxStart), sizeof(nxStart));
    outputFile.write(reinterpret_cast<char *>(&nyStart), sizeof(nyStart));
    outputFile.write(reinterpret_cast<char *>(&nzStart), sizeof(nzStart));
    outputFile.write(reinterpret_cast<char *>(&mx), sizeof(mx));
    outputFile.write(reinterpret_cast<char *>(&my), sizeof(my));
    outputFile.write(reinterpret_cast<char *>(&mz), sizeof(mz));
    outputFile.write(reinterpret_cast<char *>(&xLength), sizeof(xLength));
    outputFile.write(reinterpret_cast<char *>(&yLength), sizeof(yLength));
    outputFile.write(reinterpret_cast<char *>(&zLength), sizeof(zLength));
    outputFile.write(reinterpret_cast<char *>(&alpha), sizeof(alpha));
    outputFile.write(reinterpret_cast<char *>(&beta), sizeof(beta));
    outputFile.write(reinterpret_cast<char *>(&gamma), sizeof(gamma));
    outputFile.write(reinterpret_cast<char *>(&mapc), sizeof(mapc));
    outputFile.write(reinterpret_cast<char *>(&mapr), sizeof(mapr));
    outputFile.write(reinterpret_cast<char *>(&maps), sizeof(maps));
    outputFile.write(reinterpret_cast<char *>(&dMin), sizeof(dMin));
    outputFile.write(reinterpret_cast<char *>(&dMax), sizeof(dMax));
    outputFile.write(reinterpret_cast<char *>(&dMean), sizeof(dMean));
    outputFile.write(reinterpret_cast<char *>(&ispg), sizeof(ispg));
    outputFile.write(reinterpret_cast<char *>(&nsymbt), sizeof(nsymbt));
    outputFile.write(reinterpret_cast<char *>(extra), 25*sizeof(int));
    outputFile.write(reinterpret_cast<char *>(&xOrigin), sizeof(xOrigin));
    outputFile.write(reinterpret_cast<char *>(&yOrigin), sizeof(yOrigin));
    outputFile.write(reinterpret_cast<char *>(&zOrigin), sizeof(zOrigin));
    outputFile.write(reinterpret_cast<char *>(map), 4*sizeof(char));
    outputFile.write(reinterpret_cast<char *>(&machineStamp), sizeof(machineStamp));
    outputFile.write(reinterpret_cast<char *>(&rms), sizeof(rms));
    outputFile.write(reinterpret_cast<char *>(&nlabl), sizeof(nlabl));
    outputFile.write(reinterpret_cast<char *>(label), 10*80*sizeof(char));

    if (nx < 0 || ny < 0 || nz < 0) {
        printf("negetive nx , ny, nz value is MRC map");
        return;
    }
    int count = 0;

    for(int i = nzS; i < nzE ; i++)
    {
        for(int j = nyS; j < nyE; j++)
        {
            for(int k = nxS; k < nxE; k++)
            {
              int nxNew = ((k-c1) * cos(angle*3.14/180)) - ((j-c2) * sin(angle*3.14/180)) + c1;
              int nyNew = ((k-c1) * sin(angle*3.14/180)) + ((j-c2) * cos(angle*3.14/180)) + c2;
              if(nxNew > -1 && nxNew <= 300 && nyNew > -1 && nyNew <= 300 )
              {
                  outputFile.write(reinterpret_cast<char *>(&cube[nxNew][nyNew][i]), sizeof(float));
              }

            }
            // break;
        }
        // break;
    }

    outputFile.close();

    printf("\nwriting completed\n\n\n");
}


/*
 * Print the value of all MRC file header field
*/

void MRC::showMRCFileDescription()
{
    cout<<"nx, ny, nz: "<<nx<<","<<ny<<","<<nz<<endl;
    cout<<"mode: "<<mode<<endl;
    cout<<"nxStart, nyStart, nzStart: "<<nxStart<<","<<nyStart<<","<<nzStart<<endl;
    cout<<"mx, my, mz: "<<mx<<","<<my<<","<<mz<<endl;
    cout<<"xLength, yLength, zLength: "<<xLength<<","<<yLength<<","<<zLength<<endl;
    cout<<"alpha, beta, gamma: "<<alpha<<","<<beta<<","<<gamma<<endl;
    cout<<"mapc, mapr, maps: "<<mapc<<","<<mapr<<","<<maps<<endl;
    cout<<"dMin, dMax, dMean: "<<minDensity()<<","<<maxDensity()<<","<<meanDensity()<<endl;
    cout<<"ispg: "<<ispg<<endl;
    cout<<"nsymbt: "<<nsymbt<<endl;
    for(int i = 0; i < 25; i ++)
        cout<<extra[i]<<" ";
    cout<<endl;
    cout<<"xOrigin, yOrigin, zOrigin: "<<xOrigin<<","<<yOrigin<<","<<zOrigin<<endl;
    cout<<"map: "<<map<<endl;
    cout<<"machineStamp: "<<machineStamp<<endl;
    cout<<"rms: "<<rms<<endl;
    cout<<"nlabl: "<<nlabl<<endl;
    cout<<"-------------------------------------------"<<endl;

    cout<<endl;
}

float MRC:: minDensity()
{
    float minimumDensity = DENSITY_INFINITY;

    for(int i = 0; i < nz; i++) {
        for(int j = 0; j < ny; j++) {
            for(int k = 0; k < nx; k++) {

                if (minimumDensity > cube[k][j][i]) {

                    minimumDensity = cube[k][j][i];
                }

            }
        }
    }

    return minimumDensity;

}

float MRC:: maxDensity()
{
    float maximumDensity = -DENSITY_INFINITY;

    for(int i = 0; i < nz; i++) {
        for(int j = 0; j < ny; j++) {
            for(int k = 0; k < nx; k++) {

                if (maximumDensity < cube[k][j][i]) {

                    maximumDensity = cube[k][j][i];
                }
            }
        }
    }

    return maximumDensity;
}

float MRC:: meanDensity()
{
    float totalDensity = 0;

    for(int i = 0; i < nz; i++) {
        for(int j = 0; j < ny; j++) {
            for(int k = 0; k < nx; k++) {

                totalDensity += cube[k][j][i];
            }
        }
    }

    return (float) totalDensity / (nz * ny * nx);
}




#pragma mark -- Find Amino Acid
MRC::~MRC()
{
    for(int i = 0; i < nx; i++)
        for(int j = 0; j < ny; j++)
        {
            delete [] cube[i][j];
        }

    for(int i = 0; i < nx; i++)
        delete [] cube[i];

    delete [] cube;
}
