/******** FLAVOUR  *********//******** ACTION **********//******* MODEL *******/ 
/*  0 - light only         *  *  0 - fit               *  *   0 - gbw         * 
 *  1 - light + heavy      *  *  1 - list              *  *   1 - bgk         *
 *  2 - charm only         *  *  2 - test              *  *                   * 
 *  3 - beuaty only        *  *                        *  *                   */
/***************************//**************************//*********************/
int  flavour      = 1 ;       int  action       = 0;      int  model        = 0;
/******* DATAFORM **********//****** CHEB_APPROX *******//***** GLUON_INT *****/
/*  0 - reduced cs         *  * 0 - exact form         *  * 0 - Golec Simpson * 
 *  1 - F_2                *  * 1 - with Cheb. approx. *  * 1 - CERN Simpson  *
 *  2 - F_L                *  *                        *  *                   */
/***************************//**************************//*********************/
int  dataform     = 1;       int  cheb_approx  =1 ;      int  gluon_int    = 0;
/******* UIF_INT ***********//****** xB MODIFICATION ***//*****   BEUTY   *****/
/*  0 - dadmul             *  * 0 - without mod        *  *  0 - without b    *
 *  1 - simps2d            *  * 1 - wth mod            *  *  1 - with b       *
 *                         *  *                        *  *                   */
/***************************//**************************//*********************/
int  uif_int      = 1;       int  xbj_mod      = 1;      int  fl_beauty    = 1;
/******* DATATYPE***********//******* DATATYPE***********/
/*  0 - previous           *  *  0 - pure gbw           * 
 *  1 - new                *  *  1 - gbw with sudakov   * 
 *                         *  *                         */
/***************************//***************************/
int  datatype     = 1;       int  sudflag = 2;     


int  light_charm;

int photo = 0;

int  temp_var;
