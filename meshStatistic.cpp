/*
 * meshStatistic.cpp
 *
 * Written by Jose Miguel Espadero <josemiguel.espadero@urjc.es>
 *
 * This code is written as material for the FMF class of the
 * Master Universitario en Informatica Grafica, Juegos y Realidad Virtual.
 * Its purpose is to be didactic and easy to understand, not hard optimized.
 *
 * This file compute some statistics about a mesh and dump then to the console.
 * Also produce an output mesh with degenerate triangles remarked.
 * 
 * //TODO: Fill-in your name and email
 * Name of alumn: 
 * Email of alumn:
 * Year: 2022
 * 
 */

#ifdef _MSC_VER
#pragma warning(error: 4101)
#endif

#define _CRT_NONSTDC_NO_DEPRECATE
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <SimpleMesh.hpp>
#include <ColorMesh.hpp>

int main (int argc, char *argv[])
{
    try
    {
        //Set default input mesh filename
        std::string filename("mallas/mannequin2.ply");
        if (argc >1)
            filename = std::string(argv[1]);

        //Set default degenerate triangle shapeFactor threshold
        double shapeFactorTh = 1/(4*sqrt(3));
        if (argc > 2)
            shapeFactorTh = strtod(argv[2], nullptr);

        ///////////////////////////////////////////////////////////////////////
        //Read a mesh from given filename
        SimpleMesh mesh;
        cout << "Loading file " << filename << endl;
        mesh.readFile(filename);

        cout << "\nNum vertex: " << mesh.numVertex() << " Num triangles: " << mesh.numTriangles()
             << " Unreferenced vertex: " << mesh.checkUnreferencedVertex() << endl;


        //Compute minimum, maximum and average area per triangle, and total area for the mesh
        double minArea, maxArea, averageArea=0, totalArea;

        //Compute minimum and maximum angle for the mesh
        double minAngle, maxAngle;

        //Compute minimum, maximum and average edge length
        double minEdgeLen, maxEdgeLen, averageEdgeLen;

        //Compute minimum, maximum and average shape factor
        double minShapeFactor, maxShapeFactor, averageShapeFactor;

        //////////////////////////////////////////////////////////////////////////////////
        //TODO 1.1:
        // Compute the values for minArea, maxArea, averageArea, totalArea
        // Compute the values for minAngle, maxAngle;
        // Compute the values for minEdgeLen, maxEdgeLen, averageEdgeLen;
        // Compute the values for minShapeFactor, maxShapeFactor, averageShapeFactor;

        //END TODO 1.1

        // Dump statistics to console output
        cout << std::fixed << std::setprecision(4) <<
                "\nArea    min: " << minArea <<
                " max: " << maxArea  <<
                " average: " << averageArea  <<
                " total: " << totalArea <<
                "\nAngle   min: " << minAngle <<
                " max: " << maxAngle  <<
                "\nEdgeLen min: " << minEdgeLen <<
                " max: " << maxEdgeLen  <<
                " average: " << averageEdgeLen  <<
                "\nShapeF  min: " << minShapeFactor <<
                " max: " << maxShapeFactor  <<
                " average: " << averageShapeFactor  << endl;

        //////////////////////////////////////////////////////////////////////////////////
        //TODO OPTATIVE 1:
        //Compute Vertex area for each vertex and store in vertexAreas vector
        //Compute minimun and maximun vertex area in minVertexArea, maxVertexArea
        std::vector<double>vertexAreas;
        double sumVertexAreas = 0;
        double minVertexArea = INFINITY;
        double maxVertexArea = 0;

        //END TODO OPTATIVE 1

        //Check values for vertex areas and compute sum of areas:
        if (vertexAreas.size() != mesh.numVertex() )
        {
          cout << "VertexAreas OPTATIVE PART NOT DONE" << endl;
        }
        else
        {
          cout << "VertexA min: " << minVertexArea <<
                  " max: " << maxVertexArea  <<
                  " average: " << sumVertexAreas / mesh.numVertex()  <<
                  " total: " << sumVertexAreas << endl;
        }


        //////////////////////////////////////////////////////////////////////////////////
        //TODO 1.2:
        //Write triangles with ShapeFactor (radius / minEdge) greater than shapeFactorTh
        //Use messages formated as: "Triangle nnn has ShapeFactor xxx"

        //END TODO 1.2

        //////////////////////////////////////////////////////////////////////////////////
        //TODO 1.3:
        //Create a colorMesh where faces with ShapeFactor greater than shapeFactorTh
        //have their vertex colored in red. Save it to file named output_statistic.ply
        //and visualize it with meshlab or another external viewer.

        //END TODO 1.3

    }

    catch (const string &str) { std::cerr << "EXCEPTION: " << str << std::endl; }
    catch (const char *str) { std::cerr << "EXCEPTION: " << str << std::endl; }
    catch (std::exception& e)    { std::cerr << "EXCEPTION: " << e.what() << std::endl;  }
    catch (...) { std::cerr << "EXCEPTION (unknow)" << std::endl; }

#ifdef WIN32
    cout << "Press Return to end the program" <<endl;
    cin.get();
#else
#endif

    return 0;
}

