typedef struct paramsdata{
    int flg_parallel;
    int flg_writemat;
    int flg_showplot;
    int flg_showconv;
    int flg_gpu;
    int K;
    int Ncfg;
    double xemit;
    double yemit;
    double zemit;
    int Npupil;
    int pixelsize;
    int Mx;
    int My;
    int Mz;
    int xrange;
    int yrange;
    int numparams;
    int Nitermax;
    int fwd;
    int depth;
    int lambda;
    int lambdacentral;
    int aberrationcorrected;
    int ringradius;
    int doelevels;
    int doephasedepth;
    int readnoisestd;
    int readnoisevariance;
    int cztN;
    int cztM;
    int cztL;
    int debugmode;
    int phtrue;
    int bgtrue;
    int varfit;

    double azim;
    double m;
    double alpha;
    double tollim;
    double NA;
    double refmed;
    double refcov;
    double refimm;
    double refimmnom;
    double pola;
    double welldepth;
    double g2;

    double * aberrationsoffset;
    double phi[1000];
    double theta[1000];
    double allalpha[1000];
    double Duxstr[1000];
    double Duystr[1000];


    char fitmodel[50];
    char excitation[30];
    char ztype[50];
    char dipoletype[50];
    char doetype[50];

    int zrange[2];
    int zspread[2];
    int lambdaspread[2];
    int aberrations[9][3];
    int zonefunction[9][3];

    double _Complex **Axmt, **Bxmt, **Dxmt;
    double _Complex **Aymt, **Bymt, **Dymt;

    int ***allPSFs;

    double _Complex *** wavevector, ** wavevectorzmed, ** Waberration, ****PupilMatrix;

}paramsdata;
