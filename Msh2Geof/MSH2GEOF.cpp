/*
created on 09/09/2017
*/

#include <stdlib.h>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>
#include <iterator>
#include <sstream>
#include <cmath>
#include <limits>
const bool file_is_poly=false;
const std::string test_file_name="test.msh";

const double MAXDOUBLE=(std::numeric_limits<double>::max)();
const double MINDOUBLE=(std::numeric_limits<double>::min)();

/*
enum ELEMENT_TYPE
{
    Line_2_node=1,
    Triangle_3_node=2,
    Line_3_node=8,
    Triangle_6_node=9,
    Point_1_node=15,
    Triangle_10_node=21,
    Triangle_15_node=23,
    Triangle_21_node=25,
    Line_4_node=26,
    Line_5_node=27,
    Line_6_node=28
};
*/

struct Node_msh
{
    int n;
    double x_coord;
    double y_coord;
    double z_coord;
    Node_msh(int n, double x, double y, double z):n(n), x_coord(x), y_coord(y), z_coord(z) {}
    Node_msh() {}
};

struct Element_msh
{
    int elm_number;
    int elmtype;
    int number_of_tags;
    std::vector<int> tags;
    std::vector<int> node_list;

    int back_up;
};
Element_msh ReadElementfromMsh(std::string line);
void DetermineElementType(bool *OK2output, std::string &type_name, const int n );
bool IsEdge(int n);
std::string EdgeName(int n);

int main(int argc, char *argv[])
{
    double XMAX,XMIN,YMAX,YMIN,ZMAX,ZMIN;
    XMAX=MINDOUBLE;YMAX=MINDOUBLE;ZMAX=MINDOUBLE;
    XMIN=MAXDOUBLE;YMIN=MAXDOUBLE;ZMIN=MAXDOUBLE;
///////////////////////////////////////////////////////////////////
    std::string mshFileName,geofFileName;
    if(argc==2)
    {
        mshFileName=argv[1];
    }
    else
    {
        mshFileName=test_file_name;// editable!!!!!!!!!!!!!!!!!!!!!!!!!
    }

    mshFileName.erase(mshFileName.end()-4,mshFileName.end() );
    geofFileName=mshFileName+".geof";
    mshFileName+=".msh";

    std::ifstream ReadMSH(mshFileName.c_str(), std::ios::in);
    if(!ReadMSH)
    {
        std::cout<<"Can not open *.msh file"<<std::endl;
        exit(1);
    }
    std::cout<<"Open File Successfully"<<std::endl;

    std::vector<Node_msh> AllNodes;
    std::vector<Element_msh> AllElements;
    int number_of_nodes,number_of_elements;

    std::istringstream Insert;
    std::string linePresent;

//******************************************************************************
                ////////////////////
                // Read msh File  //
                ////////////////////
//******************************************************************************

{
//$MeshFormat
//version_number file-type data-size
//$EndMeshFormat
//$PhysicalNames
//number-of-names
//physical-dimension physical-number "physical-name"
//...
//$EndPhysicalNames
//$Nodes
    while (linePresent.find("Nodes")==std::string::npos)
    {
        std::getline( ReadMSH, linePresent);
    }
}

//number-of-nodes
{
    std::getline( ReadMSH, linePresent);
    Insert.str(linePresent);
    Insert>>number_of_nodes;
}

//node-number x-coord y-coord z-coord   //(number_of_nodes lines)
{


    for (int i=0; i< number_of_nodes; i++)
    {   int tempn;
        double tempx,tempy,tempz;
        std::getline(ReadMSH, linePresent);
        std::istringstream InsFlow;
        InsFlow.str(linePresent);
        InsFlow>>tempn >>tempx >>tempy >>tempz;
        Node_msh tempnode(tempn,tempx,tempy,tempz);
        AllNodes.push_back(tempnode);

        if(tempx<XMIN) {XMIN=tempx;}
        if(tempy<YMIN) {YMIN=tempy;}
        if(tempz<ZMIN) {ZMIN=tempz;}
        if(tempx>XMAX) {XMAX=tempx;}
        if(tempy>YMAX) {YMAX=tempy;}
        if(tempz>ZMAX) {ZMAX=tempz;}

    }
}

//$EndNodes
    std::getline( ReadMSH, linePresent);
//$Elements
    std::getline( ReadMSH, linePresent);

//number-of-elements
{
    std::getline( ReadMSH, linePresent);
    std::istringstream InsFlow;
    InsFlow.str(linePresent);
    InsFlow>>number_of_elements;
}

//elm-number elm-type number-of-tags <tags> ... node-number-list
{

    for(int i=0; i< number_of_elements; i++)
    {
        Element_msh tempelement;
        std::getline(ReadMSH,linePresent);
        // call a function to deal with different type of elements
        // Element_msh MyFunction(std::string XXX);
        tempelement=ReadElementfromMsh(linePresent);

        AllElements.push_back(tempelement);
    }


}

//EndElements (also end the reading)
    std::getline( ReadMSH, linePresent);
    if(linePresent.find("EndElements")!= std::string::npos ) std::cout<<"End Reading"<<std::endl;
    else exit(1);

//******************************************************************************
                ////////////////////
                // Write Geof File//
                ////////////////////
//******************************************************************************

    std::ofstream WriteGEOF(geofFileName.c_str(),std::ios::out);
    WriteGEOF.precision(8);

WriteGEOF<<"***geometry"<<std::endl;

//Write nodes
    WriteGEOF<<"**node"<<std::endl;
    WriteGEOF<<AllNodes.size()<<' '<<"3"<<std::endl;
    int cycletimesnode;
    int cycletimeselement;
    cycletimesnode=AllNodes.size();
    cycletimeselement=AllElements.size();
    for( int i=0; i<cycletimesnode; i++)
    {
        WriteGEOF<<AllNodes[i].n<<" "<<std::scientific<<AllNodes[i].x_coord<<' '<<AllNodes[i].y_coord<<' '<<AllNodes[i].z_coord<<std::endl;
    }

//Write elements
    WriteGEOF<<"**element"<<std::endl;

    int counter=0;
    for( int i=0; i<cycletimeselement; i++)
    {

        //determine the element type
        //only triangles can be output temporarily
        //void MyFunction(bool &OK2output, std::string typename,int type );
        bool OK2output=false;
        std::string astring;
        DetermineElementType(&OK2output,astring, AllElements[i].elmtype);

        if(OK2output != false)
        {
            counter++;
            AllElements[i].elm_number=counter;//renumber the triangle or other type
            //WriteGEOF<<counter<<' '<<astring;
            //for ( int j=0; j<AllElements[i].node_list.size();j++ ) WriteGEOF<<' '<<AllElements[i].node_list[j];
            //WriteGEOF<<std::endl;

        }

    }
    WriteGEOF<<counter<<std::endl;
    counter=0;
    for( int i=0; i<cycletimeselement; i++)
    {

        //determine the element type
        //only triangles can be output temporarily
        //void MyFunction(bool &OK2output, std::string typename,int type );
        bool OK2output=false;
        std::string astring;
        DetermineElementType(&OK2output,astring, AllElements[i].elmtype);

        if(OK2output != false)
        {
            counter++;
            //AllElements[i].elm_number=counter;//renumber the triangle or other type
            WriteGEOF<<counter<<' '<<astring;
            for ( int j=0; j<AllElements[i].node_list.size();j++ ) WriteGEOF<<' '<<AllElements[i].node_list[j];
            WriteGEOF<<std::endl;

        }

    }

//Write the groups
WriteGEOF<<"***group"<<std::endl;

// Boundary left
    WriteGEOF<<"**nset x0"<<std::endl;
    for ( int i=0; i< cycletimesnode;i++ )
    {
        if ( AllNodes[i].x_coord==XMIN ) WriteGEOF<<AllNodes[i].n<<' ';
    }
    WriteGEOF<<std::endl;
//Boundary bottom
    WriteGEOF<<"**nset y0"<<std::endl;
    for ( int i=0; i< cycletimesnode;i++ )
    {
        if ( AllNodes[i].y_coord==YMIN ) WriteGEOF<<AllNodes[i].n<<' ';
    }
    WriteGEOF<<std::endl;
//Boundary right
    WriteGEOF<<"**nset x1"<<std::endl;
    for ( int i=0; i< cycletimesnode;i++ )
    {
        if ( AllNodes[i].x_coord==XMAX ) WriteGEOF<<AllNodes[i].n<<' ';
    }
    WriteGEOF<<std::endl;

//Corner leftbottom
    WriteGEOF<<"**nset leftbottom"<<std::endl;
    for ( int i=0; i< cycletimesnode;i++ )
    {
        if ( AllNodes[i].x_coord==XMIN && AllNodes[i].y_coord==YMIN ) WriteGEOF<<AllNodes[i].n<<' ';
    }
    WriteGEOF<<std::endl;

//Corner rightbottom
    WriteGEOF<<"**nset rightbottom"<<std::endl;
    for ( int i=0; i< cycletimesnode;i++ )
    {
        if ( AllNodes[i].x_coord==XMAX && AllNodes[i].y_coord==YMIN ) WriteGEOF<<AllNodes[i].n<<' ';
    }
    WriteGEOF<<std::endl;

//Corner lefttop
    WriteGEOF<<"**nset lefttop"<<std::endl;
    for ( int i=0; i< cycletimesnode;i++ )
    {
        if ( AllNodes[i].x_coord==XMIN && AllNodes[i].y_coord==YMAX ) WriteGEOF<<AllNodes[i].n<<' ';
    }
    WriteGEOF<<std::endl;

//Corner righttop
    WriteGEOF<<"**nset righttop"<<std::endl;
    for ( int i=0; i< cycletimesnode;i++ )
    {
        if ( AllNodes[i].x_coord==XMAX && AllNodes[i].y_coord==YMAX ) WriteGEOF<<AllNodes[i].n<<' ';
    }
    WriteGEOF<<std::endl;

//faset Xone Xzero Yzero (load)
if(!file_is_poly)
{
    WriteGEOF<<"**faset Xone"<<std::endl;
    for ( int i=0; i<cycletimeselement;i++ )
    {
        int ntype=AllElements[i].elmtype;
        if ( IsEdge(ntype) != false )
        {

            int rightnode,leftnode;
            leftnode=AllElements[i].node_list[0];
            rightnode=AllElements[i].node_list[1];
            if( (AllNodes[leftnode-1].x_coord==XMAX) and (AllNodes[rightnode-1].x_coord==XMAX)   )
            //if ( IsRightEdge() !=false ) 不典型，不太好
            {

                std::string astring;
                astring=EdgeName(ntype);
                WriteGEOF<<astring;
                for( int j=0; j<AllElements[i].node_list.size(); j++) WriteGEOF<<' '<<AllElements[i].node_list[j];
                WriteGEOF<<std::endl;
            }
        }

    }

    WriteGEOF<<"**faset Xzero"<<std::endl;
    for ( int i=0; i<cycletimeselement;i++ )
    {
        int ntype=AllElements[i].elmtype;
        if ( IsEdge(ntype) != false )
        {

            int rightnode,leftnode;
            leftnode=AllElements[i].node_list[0];
            rightnode=AllElements[i].node_list[1];
            if( (AllNodes[leftnode-1].x_coord==XMIN) and (AllNodes[rightnode-1].x_coord==XMIN)   )
            //if ( IsRightEdge() !=false ) 不典型，不太好
            {

                std::string astring;
                astring=EdgeName(ntype);
                WriteGEOF<<astring;
                for( int j=0; j<AllElements[i].node_list.size(); j++) WriteGEOF<<' '<<AllElements[i].node_list[j];
                WriteGEOF<<std::endl;
            }
        }

    }

     WriteGEOF<<"**faset Yzero"<<std::endl;
    for ( int i=0; i<cycletimeselement;i++ )
    {
        int ntype=AllElements[i].elmtype;
        if ( IsEdge(ntype) != false )
        {

            int rightnode,leftnode;
            leftnode=AllElements[i].node_list[0];
            rightnode=AllElements[i].node_list[1];
            if( (AllNodes[leftnode-1].y_coord==YMIN) and (AllNodes[rightnode-1].y_coord==YMIN)   )
            //if ( IsRightEdge() !=false ) 不典型，不太好
            {

                std::string astring;
                astring=EdgeName(ntype);
                WriteGEOF<<astring;
                for( int j=0; j<AllElements[i].node_list.size(); j++) WriteGEOF<<' '<<AllElements[i].node_list[j];
                WriteGEOF<<std::endl;
            }
        }

    }
}
//？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？
// 之前旧的多晶体模型还没有适配！！！！！！因为里面没有线段的分类
//？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？

//Grains Sets output
if (file_is_poly)
{   int last_line;
    last_line=-1;
    bool not_first=false;
    for( int i=0; i<cycletimeselement; i++)
    {

        int nl=0;
        //determine the element type
        //only triangles can be output temporarily
        //void MyFunction(bool &OK2output, std::string typename,int type );
        bool OK2output=false;
        std::string astring;
        DetermineElementType(&OK2output,astring, AllElements[i].elmtype);

        if(OK2output != false)
        {
            int grain_num_phy,grain_num_ele;
            grain_num_phy=AllElements[i].tags[0];
            grain_num_ele=AllElements[i].tags[1];

            if ( grain_num_ele==grain_num_phy )
            {

                if(last_line != grain_num_ele)
                {
                    if(not_first !=false) WriteGEOF<<std::endl;
                    WriteGEOF<<"**elset grain"<<grain_num_ele<<std::endl;last_line=grain_num_ele;nl=0;
                    not_first=true;
                }
                WriteGEOF<<AllElements[i].elm_number<<' ';
                nl++;
                if(nl==20) WriteGEOF<<std::endl;

            }
        }

    }

}


WriteGEOF<<"***return"<<std::endl;
WriteGEOF.close();

std::cout<<"End Writing"<<std::endl;
return 0;
}


Element_msh ReadElementfromMsh(std::string line)
{
    Element_msh tempelement;
    std::istringstream ins;
    ins.str(line);
    int number,ntype,tagnumber,temp;
    ins>>number>>ntype>>tagnumber;
    tempelement.elm_number=number;
    tempelement.elmtype=ntype;
    tempelement.number_of_tags=tagnumber;
    for(int i=0;i<tagnumber;i++)
    {
        ins>>temp;
        tempelement.tags.push_back(temp);
    }
    while(ins>>temp)
    { tempelement.node_list.push_back(temp);}

    return tempelement;
}

void DetermineElementType(bool *OK2output, std::string &type_name, const int n )
{
    if ( (n==2) or (n==9) or (n==21) or (n==23) or (n==25)   ) *OK2output=true; //Only triangles

    if      ( (n==2) ) {type_name="c2d3r";}
    else if ( (n==9) ) {type_name="c2d6r";}
    else if ( (n==21) ) {type_name="c2d10r";}
    else if ( (n==23) ) {type_name="c2d15r";}
    else if ( (n==25) ) {type_name="c2d21r";}

}


bool IsEdge(int n)
{
    if ( (n==1) or (n==8) or (n==26) or (n==27) or (n==28)   ) return true;
    else return false;

}


std::string EdgeName(int n)
{
    std::string type_name;
    if ( (n==1) ) {type_name="l2";}
    else if ( (n==8) ) {type_name="l3";}
    else if ( (n==26) ) {type_name="l4";}
    else if ( (n==27) ) {type_name="l5";}
    else if ( (n==28) ) {type_name="l6";}
    return type_name;
}
