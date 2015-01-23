/*  
     
     Program to create an X-based calculator for 
                several flow problems
  
                    Version 2.4a

                Written by Tom Benson
              NASA Lewis Research Center

      HELP information edited by Jon Goldfarb, Purdue University

     This version uses the X-FORMS Version 075 library of GUI's 
       developed by Mark Overmars and T.C. Zhao


Notice:

THIS SOFTWARE IS PUBLIC DOMAIN. IT MAY FREELY BE COPIED AND USED IN
NON-COMMERCIAL PRODUCTS, ASSUMING PROPER CREDIT TO THE AUTHOR IS GIVEN,
BUT IT MAY NOT BE RESOLD.

THE X-FORMS LIBRARY IS COPYRIGHTED BY T.C. ZHAO AND MARK OVERMARS
WITH ALL PUBLISHED AND UNPUBLISHED RIGHTS RESERVED. PERMISSION
FOR NON-COMMERCIAL AND NON-PROFIT PURPOSES IS GRANTED. COMMERCIAL
PROGRAM DEVELOPERS CONTACT THEM FOR LICENSE ARRANGEMENT.

THIS SOFTWARE IS PROVIDED ``AS IS'' WITHOUT WARRANTY OF ANY KIND,
EITHER EXPRESSED OR IMPLIED. THE ENTIRE RISK AS TO THE QUALITY AND
PERFORMANCE OF THE SOFTWARE IS WITH YOU. SHOULD THE SOFTWARE PROVE
DEFECTIVE, YOU ASSUME THE COST OF ALL NECESSARY SERVICING, REPAIR OR
CORRECTION.

   New test -
              Add Units conversion calculator
              Changes for Xforms 075
                    introduce at_close

                                               TJB 29 Apr 96
*/

#include <stdio.h>
#include "forms.h"
#include <stdlib.h>
#include <math.h>

FL_FORM *mainform, *oform, *eform, *nform, *fform;
FL_FORM *iform, *rform, *aform, *hlpform, *diagform;
FL_OBJECT *br, *killbut, *usebut, *codebut;
FL_OBJECT *diag1, *diag2, *diag3 ;
FL_OBJECT *diag4, *diag5, *diag6 ;
FL_OBJECT *diag7, *diag8, *diag9 ;
FL_OBJECT *cho2, *ino1, *ino2, *outo, *gamo;
FL_OBJECT *che2, *ine1, *ine2, *oute, *game;
FL_OBJECT *chn1, *inn1, *outn, *gamn;
FL_OBJECT *chi1, *chi2, *ini1, *outi, *gami;
FL_OBJECT *chr1, *chr2, *inr1, *inr2, *gamr, *outr;
FL_OBJECT *chf1, *chf2, *inf1, *inf2, *gamf, *outf;
FL_OBJECT *ina1, *ina2, *outa;
float mach1,mach2;
float gama,gm1,gp1 ;
float rhorat,prat,trat,ptrat,ang,ttrat;
float alt,ps,ts,u0,q0,rey;
float delmax,thetmax,shkang,expang ;
float arat,vrat,mu,nu,p2p1,t2t1,r2r1,wcor,qrat,frat,leng ;
float delr,delmin,convdr=3.1415926/180. ;
int isub, inopt1, inopt2, erflo, gamopt;

/*       ******** FORM DEFINITIONS and UTILITIES *******************    */

/*      ------------    Utilities  ------- */
  
void fl_putr (FL_OBJECT *box, float number, int precis)
/* Utility to put floating number into output box */
{
    char str1[7];
    switch (precis) {
        case 0:  sprintf(str1,"%7.0f",number); break ;
        case 1:  sprintf(str1,"%7.1f",number); break ;
        case 2:  sprintf(str1,"%7.2f",number); break ;
        case 3:  sprintf(str1,"%7.3f",number); break ;
        case 4:  sprintf(str1,"%7.4f",number); break ;
        case 5:  sprintf(str1,"%7.5f",number); break ;
    }
    fl_set_input(box,str1);
}

float fl_getr (FL_OBJECT *box)     
/* Utility to get a floating point number from a box */
{
    float number;
    char str1[7];
    sprintf(str1,"%s",fl_get_input(box));
    number = atof (str1) ;
    return(number) ;
}

void fl_puti (FL_OBJECT *box, int number)
/* Utility to put integer into output box */
{
    char str1[5];
    sprintf(str1,"%3.1d",number);
    fl_set_input(box,str1);
}

int fl_geti (FL_OBJECT *box)
/* Utility to get an integer number from a box */
{
    int number;
    char str1[5];
    sprintf(str1,"%s",fl_get_input(box));
    number = atoi (str1) ;
    return(number) ;
}

/*      ------------    Background (MAIN) Form  ------- */

int at_close(FL_FORM *form, void *dat)
{
   return (FL_IGNORE) ;
}

void clos_swtch (FL_OBJECT *obj, long arg)
/*   switching routine to close windows */
{
  switch (arg)
  {
      case 1:  fl_hide_form(oform)    ; break ;  /* close oblique shk form */
      case 2:  fl_hide_form(eform)    ; break ;  /* close expansion form   */
      case 3:  fl_hide_form(diagform) ; break ;  /* close diag form        */
      case 4:  fl_hide_form(hlpform)  ; break ;  /* close help form        */
      case 5:  fl_hide_form(nform)    ; break ;  /* close normal shk form  */
      case 6:  fl_hide_form(iform)    ; break ;  /* close isentropic form  */
      case 7:  fl_hide_form(rform)    ; break ;  /* close rayleigh form    */
      case 8:  fl_hide_form(aform)    ; break ;  /* close atmosphere form  */
      case 9:  fl_hide_form(fform)    ; break ;  /* close fanno form    */
  }
}

void help_swtch (FL_OBJECT *obj, long arg)
/*   switching routine for help screens*/
{
  fl_clear_browser(br) ;
  switch (arg)
  {
   case 1: fl_load_browser(br,"help.gen") ; break;  /* General */
   case 2: fl_load_browser(br,"help.isen") ; break;  /* Isentropic */
   case 3: fl_load_browser(br,"help.oblq") ; break;  /* Oblique shock */
   case 4: fl_load_browser(br,"help.norm") ; break;  /* Normal Shock */
   case 5: fl_load_browser(br,"help.exp") ; break;  /* Centered Expansion */
   case 6: fl_load_browser(br,"help.ray") ; break;  /* Rayleigh Flow */
   case 7: fl_load_browser(br,"help.atm") ; break;  /* Atmosphere */
   case 8: fl_load_browser(br,"help.ads") ; break;  /* Advertisemants */
   case 9: fl_load_browser(br,"help.fan") ; break;  /* Fanno Flow */
  }
}

void create_hlpform()
{
  FL_OBJECT *obj;
  hlpform = fl_bgn_form(FL_UP_BOX,650.0,320.0);

  br = fl_add_browser(FL_NORMAL_BROWSER,10.0,10.0,500.0,300.0,"");
     fl_set_object_color(br,7,0);
  obj =fl_add_button(FL_NORMAL_BUTTON,520.0,280.0,120.0,25.0,"General");
     fl_set_call_back(obj,help_swtch,1);
  obj =fl_add_button(FL_NORMAL_BUTTON,520.0,245.0,120.0,25.0,"Isentropic");
     fl_set_call_back(obj,help_swtch,2);
  obj =fl_add_button(FL_NORMAL_BUTTON,520.0,215.0,120.0,25.0,"Oblique Shock");
     fl_set_call_back(obj,help_swtch,3);
  obj =fl_add_button(FL_NORMAL_BUTTON,520.0,185.0,120.0,25.0,"Normal Shock");
     fl_set_call_back(obj,help_swtch,4);
  obj =fl_add_button(FL_NORMAL_BUTTON,520.0,155.0,120.0,25.0,"Centered Expansion");
     fl_set_call_back(obj,help_swtch,5);
  obj =fl_add_button(FL_NORMAL_BUTTON,520.0,125.0,120.0,25.0,"Rayleigh");
     fl_set_call_back(obj,help_swtch,6);
  obj =fl_add_button(FL_NORMAL_BUTTON,520.0,95.0,120.0,25.0,"Fanno");
     fl_set_call_back(obj,help_swtch,9);
  obj =fl_add_button(FL_NORMAL_BUTTON,520.0,65.0,120.0,25.0,"Atmosphere");
     fl_set_call_back(obj,help_swtch,7);
  obj =fl_add_button(FL_NORMAL_BUTTON,520.0,35.0,120.0,25.0,"Ads");
     fl_set_call_back(obj,help_swtch,8);
  obj =fl_add_button(FL_NORMAL_BUTTON,520.0,5.0,120.0,25.0,"Close");
     fl_set_object_color(obj,3,3);
     fl_set_call_back(obj,clos_swtch,4);
  fl_end_form();
}

void main_swtch (FL_OBJECT *obj, long arg) 
/*   switching routine for main - uses choice button*/
{
  int i1 ;

  i1 = fl_get_menu(codebut) ;

  switch (i1)
  {
      case 1: {                                /* Isentropic */
            fl_show_form(iform,FL_PLACE_ASPECT,FL_FULLBORDER,"Isentropic");
            break ;
      }
      case 2: {                                  /* Normal Shock */
            fl_show_form(nform,FL_PLACE_ASPECT,FL_FULLBORDER,"N-Shock");
            break ; 
      }
      case 3: {                                 /* Oblique Shock */
            fl_show_form(oform,FL_PLACE_ASPECT,FL_FULLBORDER,"Oblique Shock");
            break;
      }
      case 4: {                                  /* Expansion Fan */
            fl_show_form(eform,FL_PLACE_ASPECT,FL_FULLBORDER,"Expansion");
            break ;
      }
      case 5: {                                /* Rayleigh Flow */
            fl_show_form(rform,FL_PLACE_ASPECT,FL_FULLBORDER,"Rayleigh Flow");
            break ;
      }
      case 6: {                                /* Fanno Flow */
            fl_show_form(fform,FL_PLACE_ASPECT,FL_FULLBORDER,"Fanno Flow");
            break ;
      }
      case 7: {                                /* Atmospheric  */
            fl_show_form(aform,FL_PLACE_ASPECT,FL_FULLBORDER,"Atmos");
            break ;
      }
  }
}  

void create_diagform()
{
  FL_OBJECT *obj;
  diagform = fl_bgn_form(FL_UP_BOX,350.0,120.0);

  diag1 = obj =fl_add_input(FL_NORMAL_INPUT,55.,90.,60.,30.,"diag1");
         fl_set_object_color(obj,7,7);
  diag2 = obj =fl_add_input(FL_NORMAL_INPUT,155.,90.,60.,30.,"diag2");
         fl_set_object_color(obj,7,7);
  diag3 = obj =fl_add_input(FL_NORMAL_INPUT,255.,90.,60.,30.,"diag3");
         fl_set_object_color(obj,7,7);
  diag4 = obj =fl_add_input(FL_NORMAL_INPUT,55.,60.,60.,30.,"diag4");
         fl_set_object_color(obj,7,7);
  diag5 = obj =fl_add_input(FL_NORMAL_INPUT,155.,60.,60.,30.,"diag5");
         fl_set_object_color(obj,7,7);
  diag6 = obj =fl_add_input(FL_NORMAL_INPUT,255.,60.,60.,30.,"diag6");
         fl_set_object_color(obj,7,7);
  diag7 = obj =fl_add_input(FL_NORMAL_INPUT,55.,30.,60.,30.,"diag7");
         fl_set_object_color(obj,7,7);
  diag8 = obj =fl_add_input(FL_NORMAL_INPUT,155.,30.,60.,30.,"diag8");
         fl_set_object_color(obj,7,7);
  diag9 = obj =fl_add_input(FL_NORMAL_INPUT,255.,30.,60.,30.,"diag9");
         fl_set_object_color(obj,7,7);
  obj =fl_add_button(FL_NORMAL_BUTTON,130.0,5.0,100.0,25.0,"Close");
     fl_set_object_color(obj,3,3);
     fl_set_call_back(obj,clos_swtch,3);
  fl_end_form();
}

void usept_swtch (FL_OBJECT *obj, long arg)
/*   switching routine for main screen buttons */
{
  int i1 ;

  i1 = fl_get_menu(usebut) ;
  switch (i1)
  {
   case 1:
      {
         fl_show_form(hlpform,FL_PLACE_ASPECT,FL_FULLBORDER,"Help");
         fl_load_browser(br,"help.one") ;
         break ;
      }
   case 2:
      {
         fl_show_form(diagform,FL_PLACE_ASPECT,FL_FULLBORDER,"Diagnostic");
         break ;
      }
  }
}

void create_mainform()
{
  FL_OBJECT *obj;

  mainform = fl_bgn_form(FL_UP_BOX,200.0,260.0);
  obj = fl_add_box(FL_FLAT_BOX,5.0,5.0,190.0,250.0,"");
    fl_set_object_color(obj,7,7);

               /* Label */
    obj = fl_add_box(FL_NO_BOX,20.0,215.0,150.0,30.0,"Flow Calculators");
      fl_set_object_lsize(obj,FL_LARGE_SIZE);
      fl_set_object_lstyle(obj,1);
      fl_set_object_lcol(obj,4);
    obj = fl_add_box(FL_NO_BOX,20.0,180.0,150.0,25.0,"NASA Lewis");
      fl_set_object_lsize(obj,FL_MEDIUM_SIZE);
      fl_set_object_lstyle(obj,1);
      fl_set_object_lcol(obj,1);
    obj = fl_add_box(FL_NO_BOX,20.0,145.0,150.0,25.0,"Version 2.4");
      fl_set_object_lstyle(obj,1);
      fl_set_object_lsize(obj,FL_NORMAL_SIZE);

  usebut=obj=fl_add_menu(FL_PUSH_MENU,15.0,90.0,170.0,30.0,"User Options");
                fl_set_object_lsize(obj,FL_MEDIUM_SIZE);
                fl_set_object_boxtype(obj,FL_RSHADOW_BOX);
                fl_addto_menu(obj,"Help");
                fl_addto_menu(obj,"Diagnostics");
                fl_set_call_back(obj,usept_swtch,0);

  codebut=obj=fl_add_menu(FL_PUSH_MENU,15.0,50.0,170.0,30.0,"Selections");
                fl_set_object_lsize(obj,FL_MEDIUM_SIZE);
                fl_set_object_boxtype(obj,FL_RSHADOW_BOX);
                fl_addto_menu(obj,"Isentropic Flow");
                fl_addto_menu(obj,"Normal Shock");
                fl_addto_menu(obj,"Oblique Shock");
                fl_addto_menu(obj,"Centered Expansion");
                fl_addto_menu(obj,"Rayleigh Flow");
                fl_addto_menu(obj,"Fanno Flow");
                fl_addto_menu(obj,"Standard Atmosphere");
                fl_set_call_back(obj,main_swtch,0);

              /* Exit Button */
  killbut = obj = fl_add_button(FL_NORMAL_BUTTON,55.0,5.0,80.0,25.0,"Exit");
    fl_set_object_color(obj,1,1);
    fl_set_object_lcol(obj,7);
  fl_end_form();
}

/*      ------------    Isentropic Flow Calculator   ------- */

void ise_swtch (FL_OBJECT *obj, long arg) 
/*   callback routine for Isentropic Flow */
{
  void anlyi () ;

  gama = fl_getr(gami) ;          /* get gamma */
  anlyi () ;
  if (erflo != 1) {
    switch (arg)
    {                   /* output */
     case 1: fl_putr (outi,mach1,4) ; break ;  /*  Mach */
     case 2: fl_putr (outi,wcor,4); break ;    /* Corrected Airflow/Area */
     case 3: fl_putr (outi,qrat,4); break ;    /* q / p */
     case 4: fl_putr (outi,prat,4) ; break ;   /* pressure ratio */
     case 5: fl_putr (outi,arat,4) ;  break ;  /* area ratio */
     case 6: fl_putr (outi,rhorat,4); break ;  /* density */
     case 7: fl_putr (outi,trat,4) ;   break ; /* temperature */
     case 8: fl_putr (outi,vrat,4) ; break ;   /* v over a star */
     case 9: fl_putr (outi,nu,4) ;  break ;    /* P-M function */
     case 10: fl_putr (outi,mu,4); break ;     /* Mach Angle */
    }
  }
}   

void create_iform()
{
  FL_OBJECT *obj;

  iform = fl_bgn_form(FL_UP_BOX,200.0,280.0);

               /* Input */
    chi1 = fl_add_choice(FL_NORMAL_CHOICE,5.0,250.0,90.0,25.0,"");
       fl_set_object_boxtype(chi1,FL_RSHADOW_BOX);
       fl_addto_choice(chi1,"Mach ");
       fl_addto_choice(chi1,"Wc/A");
       fl_addto_choice(chi1,"q/p");
       fl_addto_choice(chi1,"p/pt");
       fl_addto_choice(chi1,"A/A*");
       fl_addto_choice(chi1,"r/rt");
       fl_addto_choice(chi1,"T/Tt");
       fl_addto_choice(chi1,"V/a*");
       fl_addto_choice(chi1,"P-M");
       fl_addto_choice(chi1,"M ang");
    ini1 = fl_add_input(FL_NORMAL_INPUT,5.,220.,90.,25.,"") ;
       fl_set_object_color(ini1,7,7);
       fl_set_object_lsize(ini1,FL_MEDIUM_SIZE); 
             /* gamma window */
    obj = fl_add_box(FL_NO_BOX,100.0,250.0,90.0,25.0,"g");
       fl_set_object_lstyle(obj,15) ;
       fl_set_object_lsize(obj,FL_NORMAL_SIZE);
    gami = obj = fl_add_input(FL_NORMAL_INPUT,100.,220.,90.,25.,"") ;
       fl_set_object_color(obj,7,7); 
       fl_set_object_lsize(obj,FL_MEDIUM_SIZE); 


             /* Choice of Outputs */
    obj = fl_add_box(FL_FLAT_BOX,5.0,30.0,190.0,185.0,"");
             fl_set_object_color(obj,7,7);
       obj = fl_add_lightbutton(FL_RADIO_BUTTON,5.0,185.0,90.0,25.0,"Mach");
         fl_set_object_callback(obj,ise_swtch,1); 
       obj = fl_add_lightbutton(FL_RADIO_BUTTON,100.0,185.0,90.0,25.0,"Wc/A");
         fl_set_object_callback(obj,ise_swtch,2); 
       obj = fl_add_lightbutton(FL_RADIO_BUTTON,5.0,155.0,90.0,25.0,"q/p");
         fl_set_object_callback(obj,ise_swtch,3); 
       obj = fl_add_lightbutton(FL_RADIO_BUTTON,100.0,155.0,90.0,25.0,"p/pt");
         fl_set_object_callback(obj,ise_swtch,4); 
       obj = fl_add_lightbutton(FL_RADIO_BUTTON,100.0,125.0,90.0,25.0,"A/A*");
         fl_set_object_callback(obj,ise_swtch,5); 
       obj = fl_add_lightbutton(FL_RADIO_BUTTON,5.0,125.0,90.0,25.0,"r/rt");
         fl_set_object_callback(obj,ise_swtch,6); 
       obj = fl_add_lightbutton(FL_RADIO_BUTTON,5.0,95.0,90.0,25.0,"T/Tt");
         fl_set_object_callback(obj,ise_swtch,7); 
       obj = fl_add_lightbutton(FL_RADIO_BUTTON,100.0,95.0,90.0,25.0,"V/a*");
         fl_set_object_callback(obj,ise_swtch,8); 
       obj = fl_add_lightbutton(FL_RADIO_BUTTON,5.0,65.0,90.0,25.0,"P-M ");
         fl_set_object_callback(obj,ise_swtch,9); 
       obj = fl_add_lightbutton(FL_RADIO_BUTTON,100.0,65.0,90.0,25.0,"M ang");
         fl_set_object_callback(obj,ise_swtch,10); 
    chi2 = fl_add_choice(FL_NORMAL_CHOICE,7.0,33.0,80.0,25.0,"");
       fl_set_object_boxtype(chi2,FL_RSHADOW_BOX);
       fl_addto_choice(chi2,"Super");
       fl_addto_choice(chi2,"Sub");
             /* Output window */
   outi = fl_add_input(FL_NORMAL_INPUT,95.,35.,100.,25.,"") ;
       fl_set_input_color(outi,0,7);
       fl_set_object_lsize(outi,FL_MEDIUM_SIZE);

              /* Exit Button */
    obj = fl_add_button(FL_NORMAL_BUTTON,55.0,5.0,80.0,25.0,"Exit");
         fl_set_object_callback(obj,clos_swtch,6); 
         fl_set_object_color(obj,1,1);
         fl_set_object_lcol(obj,7);

  fl_end_form();
}

void anlyi ()                           /* analysis for isentropic flow */
{
   float fac;         
   float get_mach(float, int);
   void get_machpm () ;
   void get_mach_arat ();
   void isen () ;
   
   erflo = 0 ;
   if (gama <1.01) gama = 1.01 ;
   gm1 = gama - 1.0 ;
   gp1 = gama + 1.0;
   inopt1 = fl_get_choice(chi1);
   isub = fl_get_choice(chi2) ;

    switch (inopt1) {
        case 1: {                               /* mach number given */
             mach1 = fl_getr(ini1) ;
             break ;
        }
        case 2: {                               /* Corrected Airflow given */
             wcor = fl_getr(ini1) ;
             mach1 = get_mach(wcor,isub) ;
             break ;
        }
        case 3: {                               /* q ratio given */
             qrat = fl_getr(ini1) ;
             mach1 = sqrt(2.0*qrat/gama) ;
             break ;
        }
        case 4: {                               /* pressure ratio given */
             prat = fl_getr(ini1) ;
             mach1 = sqrt(2.0*(pow(prat,-gm1/gama) - 1.0)/gm1) ;
             break ;
         }
        case 5: {                               /* area ratio given */
             arat = fl_getr(ini1) ;
             get_mach_arat() ;
             break ;
         }
        case 6: {                               /* density ratio given */
             rhorat = fl_getr(ini1) ;
             mach1 = sqrt(2.0*(pow(rhorat,-gm1) - 1.0)/gm1) ;
             break ;
         }
        case 7: {                               /* temperature ratio given */
             trat = fl_getr(ini1) ;
             mach1 = sqrt(2.0*(pow(trat,-1.0) - 1.0)/gm1) ;
             break ;
         }
        case 8: {                               /* velocity ratio given */
             vrat = fl_getr(ini1) ;
             mach1 = sqrt(2.0*vrat*vrat/(gp1*(1.0-gm1*vrat*vrat/gp1))) ;
             break ;
         }
        case 9: {                               /* Prandtl-Meyer given */
             nu = fl_getr(ini1) ;
             get_machpm() ;
             break ;
         }
        case 10: {                               /* Mach angle given */
             mu = fl_getr(ini1) ;
             if (isub != 1) {erflo =1; fl_set_input(outi,"Subsonic"); return ;}
             else {mach1 = 1.0/sin(mu*convdr) ; break; }
         }
     }
     isen ();
     return ;
}

float get_air(float mach)
/* Utility to get the corrected airflow per area given the Mach number */
{
    float number, fac1;
    fac1 = pow((1.0+.2*mach*mach),3.0);
    number =  .59352 * mach/ fac1 ;

    return(number) ;
}

float get_mach(float corair, int super)
/* Utility to get the Mach number given the corrected airflow per area */
{
    float number;                     /* iterate for mach number */
    float fac1,deriv,machn,macho,airo,airn;

    if (corair > .3432) {
      return ;
    }

      airo = .25618 ;                 /* initial guess */
      if (super == 1) macho = 1.703 ; /* supersonic */
      else macho = .5;                /* subsonic */
      machn = macho - .2  ;
      while (fabs(corair - airo) > .0001) {
         airn =  get_air(machn) ;
         deriv = (airn-airo)/(machn-macho) ;
         airo = airn ;
         macho = machn ;
         machn = macho + (corair - airo)/deriv ;
      }
      number = macho ;

    return(number) ;
}

void get_mach_arat()                       /*  get the Mach number */
{                                          /* given the area ratio A/Astar */
    float deriv,machn,macho,aro,arn,fac,fac1;

      fac1 = gp1/(2.0*gm1) ;

      aro = 2.0 ;                 /* initial guess */
      macho = 2.2 ;                      /* supersonic */
      if (isub != 1) macho = .30;        /* subsonic */
      machn = macho + .05  ;
      while (fabs(arat - aro) > .0001) {
         fac = 1.0 + .5*gm1*machn*machn;
         arn = 1.0/(machn * pow(fac,-fac1) * pow(gp1/2.0,fac1))  ;
         deriv = (arn-aro)/(machn-macho) ;
         aro = arn ;
         macho = machn ;
         machn = macho + (arat - aro)/deriv ;
      }
      mach1 = macho ;
      return ;
}

float get_pisen(float machin, float gama)
/* Utility to get the isentropic pressure ratio given the mach number */
{
    float number,fac1,gm1,mach1s;
    mach1s = machin*machin ;
    gm1 = gama - 1.0 ;
    fac1 = 1.0 + .5*gm1*mach1s;
    number = pow(1.0/fac1,gama/gm1) ;

    return(number) ;
}

float get_tisen(float machin, float gama)
/* Utility to get the isentropic temperature ratio given the mach number */
{
    float number,gm1,mach1s;
    mach1s = machin*machin ;
    gm1 = gama - 1.0 ;
    number = 1.0 /(1.0 + .5*gm1*mach1s) ;

    return(number) ;
}

void get_machpm ()                        /* get the Mach number */
{                                         /* given the Prandtl-Meyer angle */
   float msm1o,msm1n;
   float nur,nuo,nun,deriv ;

   nur = nu*convdr ;

   msm1o = 2.0;                                  /* iterate for mach */
   nuo = sqrt(gp1/gm1)*atan(sqrt(gm1*msm1o/gp1)) - atan(sqrt(msm1o));
   msm1n = msm1o+.1 ;
   while (fabs(nur - nuo) > .0001) {
         nun = sqrt(gp1/gm1)*atan(sqrt(gm1*msm1n/gp1)) - atan(sqrt(msm1n));
         deriv = (nun-nuo)/(msm1n-msm1o) ;
         nuo = nun ;
         msm1o = msm1n ;
         msm1n = msm1o + (nur-nuo)/deriv ;
    }
   mach1 = sqrt(msm1o + 1.0);
   return ;
}

void isen ()                                     /* isentropic relations */  
{                                                /* from NACA 1135       */
   float mach1s,msm1,fac,fac1 ;          /* mach is given        */
                                                 /* prat and trat are ratiod */
   mach1s = mach1*mach1 ;                        /* to total conditions - */
   msm1 = mach1s - 1.0;
   fac = 1.0 + .5*gm1*mach1s;
  
   prat = pow(1.0/fac,gama/gm1) ;                       /* EQ 44 */
   trat = 1.0 / fac ;                                   /* EQ 43 */
   rhorat = pow(1.0/fac,1.0/gm1) ;                      /* EQ 45 */
   mu = (asin(1.0/mach1))/convdr ;
   nu = sqrt(gp1/gm1)*atan(sqrt(gm1*msm1/gp1)) - atan(sqrt(msm1)) ;
   nu = nu/convdr ;
   vrat = sqrt(gp1 * mach1s * trat / 2.0) ;             /* EQ 50 */
   fac1 = gp1/(2.0*gm1) ;
   arat = mach1 * pow(fac,-fac1) * pow(gp1/2.0,fac1) ;  /* EQ 80 */
   arat = 1.0/arat ;
   qrat = gama*mach1s*.5 ;                              /* EQ 47 */
   wcor = get_air(mach1) ;
   return;
}

/*      ------------    Normal Shock  Calculator  ------- */

void nor_swtch (FL_OBJECT *obj, long arg) 
/*   callback routine for Normal Shock */
{
  void anlyn () ;

  gama = fl_getr(gamn) ;          /* get gamma */
  anlyn () ;
  if (erflo !=1 ) {
     switch (arg)
     {                /* output */
      case 1: fl_putr (outn,mach1,4) ; break ;  /* upstream Mach */
      case 2: fl_putr (outn,mach2,4) ; break ;  /* downstream mach */
      case 3: fl_putr (outn,prat,4) ;  break ;  /* static pressure */
      case 4: fl_putr (outn,ptrat,4) ; break ;  /* total pressure */
      case 5: fl_putr (outn,trat,4) ;  break ;  /* temperature */
      case 6: fl_putr (outn,rhorat,4); break ;  /* density */
     }
  }
}   

void create_nform()
{
  FL_OBJECT *obj;

  nform = fl_bgn_form(FL_UP_BOX,200.0,220.0);

               /* Input */
  chn1 = fl_add_choice(FL_NORMAL_CHOICE,5.0,190.0,90.0,25.0,"");
     fl_set_object_boxtype(chn1,FL_RSHADOW_BOX);
     fl_addto_choice(chn1,"Mach 1");
     fl_addto_choice(chn1,"Mach 2");
     fl_addto_choice(chn1,"pt2/pt1");
     fl_addto_choice(chn1,"p2/p1");
     fl_addto_choice(chn1,"T2/T1");
     fl_addto_choice(chn1,"r2/r1");
  inn1 = fl_add_input(FL_NORMAL_INPUT,5.,160.,90.,25.,"") ;
     fl_set_object_color(inn1,7,7);
     fl_set_object_lsize(inn1,FL_MEDIUM_SIZE);
             /* gamma window */
    obj = fl_add_box(FL_NO_BOX,100.0,190.0,90.0,25.0,"g");
       fl_set_object_lstyle(obj,15) ;
       fl_set_object_lsize(obj,FL_NORMAL_SIZE);
    gamn = obj = fl_add_input(FL_NORMAL_INPUT,100.,160.,90.,25.,"") ;
       fl_set_object_color(obj,7,7);
       fl_set_object_lsize(obj,FL_MEDIUM_SIZE);

             /* Choice of Outputs */
  obj = fl_add_box(FL_FLAT_BOX,5.0,30.0,190.0,125.0,"");
           fl_set_object_color(obj,7,7);
     obj = fl_add_lightbutton(FL_RADIO_BUTTON,5.0,125.0,90.0,25.0,"Mach 1");
       fl_set_object_callback(obj,nor_swtch,1); 
     obj = fl_add_lightbutton(FL_RADIO_BUTTON,100.0,125.0,90.0,25.0,"Mach 2");
       fl_set_object_callback(obj,nor_swtch,2); 
     obj = fl_add_lightbutton(FL_RADIO_BUTTON,5.0,95.0,90.0,25.0,"p2/p1");
       fl_set_object_callback(obj,nor_swtch,3); 
     obj = fl_add_lightbutton(FL_RADIO_BUTTON,100.0,95.0,90.0,25.0,"pt2/pt1");
       fl_set_object_callback(obj,nor_swtch,4); 
     obj = fl_add_lightbutton(FL_RADIO_BUTTON,5.0,65.0,90.0,25.0,"T2/T1");
       fl_set_object_callback(obj,nor_swtch,5); 
     obj = fl_add_lightbutton(FL_RADIO_BUTTON,100.0,65.0,90.0,25.0,"r2/r1");
       fl_set_object_callback(obj,nor_swtch,6); 
             /* Output window */
  outn = fl_add_input(FL_NORMAL_INPUT,55.,35.,100.,25.,"") ;
       fl_set_input_color(outn,0,7);
       fl_set_object_lsize(outn,FL_MEDIUM_SIZE);

              /* Exit Button */
  obj = fl_add_button(FL_NORMAL_BUTTON,55.0,5.0,80.0,25.0,"Exit");
         fl_set_object_callback(obj,clos_swtch,5); 
         fl_set_object_color(obj,1,1);
         fl_set_object_lcol(obj,7);
  fl_end_form();
}

void anlyn ()                           /* analysis for normal shock */
{
    void get_mach_ptrat () ;
    void get_mach_trat () ;
    void norms () ;

    erflo = 0;
    if (gama <1.01) gama = 1.01 ;
    gm1 = gama - 1.0 ;
    gp1 = gama + 1.0;
    inopt1 = fl_get_choice(chn1);
    shkang = 90.0 ;
    switch (inopt1) {
        case 1: {                                /* mach number given */
             mach1 = fl_getr(inn1) ;
             if (mach1 <= 1.0) {erflo=1; fl_set_input(outn,"M < 1"); return; }
             break ;
        }
        case 2: {                          /* downstream mach number given */ 
             mach2 = fl_getr(inn1) ;
             if (mach2 >= 1.0) {erflo=1; fl_set_input(outn,"M > 1"); return; }
             mach1 = sqrt((2.0+gm1*mach2*mach2)/(2.0*gama*mach2*mach2-gm1)) ;
             break ;
         }
        case 3: {                             /* total press ratio given */
             ptrat = fl_getr(inn1) ;
             if (ptrat>=1.0){erflo=1; fl_set_input(outn,"pt2/pt1 > 1"); return;}
             get_mach_ptrat() ;
             break ;
         }
        case 4: {                             /* static pressure ratio given */
             prat = fl_getr(inn1) ;
             if (prat <= 1.0) {erflo=1; fl_set_input(outn,"p2/p1 < 1"); return;}
             mach1 = sqrt((gp1 * prat + gm1)/(2.0*gama)) ;
             break ;
         }
        case 5: {                               /* temperature ratio given */
             trat = fl_getr(inn1) ;
             if (trat <= 1.0) {erflo=1; fl_set_input(outn,"T2/T1 < 1"); return;}
             get_mach_trat() ;
             break ;
         }
        case 6: {                               /* density ratio given */
             rhorat = fl_getr(inn1) ;
             if (rhorat<= 1.0){erflo=1; fl_set_input(outn,"r2/r1 < 1"); return;}
             mach1 = sqrt(2.0*rhorat/(gp1 - gm1*rhorat)) ;
             break ;
         }
     }
     norms ();
     return ;
}

void get_mach_ptrat ()                        /* get the Mach number */
{                                     /* given the total pressure ratio */
   float mso,msn;                     /* across a normal shock          */
   float ptro,ptrn,deriv ;

   mso  = 2.25;                                  /* iterate for mach */
   ptro = .9298;
   msn = mso +.2 ;
   while (fabs(ptrat - ptro) > .0001) {
         ptrn =(pow(gp1*msn/(gm1*msn + 2.0),(gama/gm1)))*
                     (pow(gp1/(2.0*gama*msn-gm1),(1.0/gm1)));
         deriv = (ptrn-ptro)/(msn-mso) ;
         ptro = ptrn ;
         mso = msn ;
         msn = mso + (ptrat-ptro)/deriv ;
    }
   mach1 = sqrt(mso);
   return ;
}

void get_mach_trat ()                        /* get the Mach number */
{                                     /* given the temperature ratio */
   float mso,msn;                     /* across a normal shock       */
   float tro,trn,deriv ;

   mso  = 2.25;                                  /* iterate for mach */
   tro = 1.320;
   msn = mso+.1 ;
   while (fabs(trat - tro) > .0001) {
         trn =((2.0*gama * msn - gm1) / gp1) /
                     (gp1 * msn / (gm1 * msn + 2.0));
         deriv = (trn-tro)/(msn-mso) ;
         tro = trn ;
         mso = msn ;
         msn = mso + (trat-tro)/deriv ;
    }
   mach1 = sqrt(mso);
   return ;
}

void norms ()                                /* normal shock relations   */
{                                            /* mach 1 given             */
   float mach2s,mach1s ;                     /* equations from NACA 1135 */

   if (mach1 < 1.0) {                  /* check for supersonic conditions */
      return ;
   }

   mach1s = mach1 * mach1 ;

   prat = (2.0 * gama * mach1s - gm1) / gp1 ;            /* EQ. 93 */
   rhorat = gp1 * mach1s / (gm1 * mach1s + 2.0) ;        /* EQ. 94 */
   trat = prat / rhorat ;                                /* EQ. 95 */
   mach2s = mach1s / (prat * rhorat) ;                   /* EQ. 96 */
   mach2 = sqrt (mach2s) ;
   ptrat =(pow(rhorat,(gama/gm1)))*
                     (pow((1.0/prat),(1.0/gm1)));       /* EQ. 99 */

   return;
}

/*      ------------    Oblique Shock  Calculator  ------- */

void oblq_swtch (FL_OBJECT *obj, long arg)
/*   callback routine for oblique shock */
{
  void anlyo () ;

  gama = fl_getr(gamo) ;          /* get gamma */
  anlyo () ;
  if (erflo !=1) {
     switch (arg)
     {                /* output */
      case 1: fl_putr (outo,mach1,4) ; break ;  /* upstream Mach */
      case 2: fl_putr (outo,mach2,4) ; break ;  /* downstream mach */
      case 3: fl_putr (outo,ang,2) ;   break ;  /* ramp angle */
      case 4: fl_putr (outo,shkang,2); break ;  /* shock angle */
      case 5: fl_putr (outo,prat,4) ;  break ;  /* static pressure */
      case 6: fl_putr (outo,ptrat,4) ; break ;  /* total pressure */
      case 7: fl_putr (outo,trat,4) ;  break ;  /* temperature */
      case 8: fl_putr (outo,rhorat,4); break ;  /* density */
      case 9: fl_putr (outo,delmax,4); break ;  /* del max */
     }
  }
}

void create_oform()
{
  FL_OBJECT *obj;

  oform = fl_bgn_form(FL_UP_BOX,290.0,220.0);

               /* Input */
  obj = fl_add_box(FL_NO_BOX,5.0,190.0,90.0,25.0,"Mach 1");
  ino1 = fl_add_input(FL_NORMAL_INPUT,5.,160.,90.,25.,"") ;
     fl_set_object_color(ino1,7,7);
     fl_set_object_lsize(ino1,FL_MEDIUM_SIZE);
  cho2 = fl_add_choice(FL_NORMAL_CHOICE,100.0,190.0,85.0,25.0,"");
     fl_set_object_boxtype(cho2,FL_RSHADOW_BOX);
     fl_addto_choice(cho2,"Ramp");
     fl_addto_choice(cho2,"Shock");
     fl_addto_choice(cho2,"pt2/pt1");
     fl_addto_choice(cho2,"p2/p1");
     fl_addto_choice(cho2,"T2/T1");
     fl_addto_choice(cho2,"r2/r1");
     fl_addto_choice(cho2,"Mach 2");
  ino2 = fl_add_input(FL_NORMAL_INPUT,100.,160.,85.,25.,"") ;
     fl_set_object_color(ino2,7,7);
     fl_set_object_lsize(ino2,FL_MEDIUM_SIZE);
             /* gamma window */
    obj = fl_add_box(FL_NO_BOX,195.0,190.0,90.0,25.0,"g");
       fl_set_object_lstyle(obj,15) ;
       fl_set_object_lsize(obj,FL_NORMAL_SIZE);
    gamo = obj = fl_add_input(FL_NORMAL_INPUT,195.,160.,90.,25.,"") ;
       fl_set_object_color(obj,7,7);
       fl_set_object_lsize(obj,FL_MEDIUM_SIZE);

             /* Choice of Outputs */
  obj = fl_add_box(FL_FLAT_BOX,5.0,30.0,285.0,125.0,"");
    fl_set_object_color(obj,7,7);
     obj = fl_add_lightbutton(FL_RADIO_BUTTON,5.0,125.0,90.0,25.0,"Mach 1");
       fl_set_object_callback(obj,oblq_swtch,1);
     obj = fl_add_lightbutton(FL_RADIO_BUTTON,100.0,125.0,90.0,25.0,"Mach 2");
       fl_set_object_callback(obj,oblq_swtch,2);
     obj = fl_add_lightbutton(FL_RADIO_BUTTON,195.0,125.0,90.0,25.0,"del max");
       fl_set_object_callback(obj,oblq_swtch,9); 
     obj = fl_add_lightbutton(FL_RADIO_BUTTON,5.0,95.0,90.0,25.0,"Ramp");
       fl_set_object_callback(obj,oblq_swtch,3);
     obj = fl_add_lightbutton(FL_RADIO_BUTTON,100.0,95.0,90.0,25.0,"Shock");
       fl_set_object_callback(obj,oblq_swtch,4);
     obj = fl_add_lightbutton(FL_RADIO_BUTTON,195.0,95.0,90.0,25.0,"p2/p1");
       fl_set_object_callback(obj,oblq_swtch,5);
     obj = fl_add_lightbutton(FL_RADIO_BUTTON,5.0,65.0,90.0,25.0,"pt2/pt1");
       fl_set_object_callback(obj,oblq_swtch,6);
     obj = fl_add_lightbutton(FL_RADIO_BUTTON,100.0,65.0,90.0,25.0,"T2/T1");
       fl_set_object_callback(obj,oblq_swtch,7);
     obj = fl_add_lightbutton(FL_RADIO_BUTTON,195.0,65.0,90.0,25.0,"r2/r1");
       fl_set_object_callback(obj,oblq_swtch,8);
             /* Output window */
  outo = fl_add_input(FL_NORMAL_INPUT,100.,35.,100.,25.,"") ;
       fl_set_input_color(outo,0,7);
       fl_set_object_lsize(outo,FL_MEDIUM_SIZE);

              /* Exit Button */
  obj = fl_add_button(FL_NORMAL_BUTTON,100.0,5.0,90.0,25.0,"Exit");
         fl_set_object_callback(obj,clos_swtch,1); 
         fl_set_object_color(obj,1,1);
         fl_set_object_lcol(obj,7);
  fl_end_form();
}

void anlyo ()                           /* analysis for oblique shock */
{
   float adet,pdet,tdet,rdet,ptdet,m2det;
   float pmax,tmax,rmax,ptmin,m2min,msave,msint;
    void get_shang_ang(float) ;
    void get_shang_mach2() ;
    void get_mach_ptrat () ;
    void get_mach_trat () ;
    void anglim() ;  
    void compr () ;

    erflo = 0 ;
    if (gama <1.01) gama = 1.01 ;
    gm1 = gama - 1.0 ;
    gp1 = gama + 1.0;
    inopt2 = fl_get_choice(cho2);
    mach1 = fl_getr(ino1);          /* get mach number*/
    if (mach1 < 1.0) {erflo=1; fl_set_input(outo,"M1 < 1"); return;}

/*      find the detachment shock angle  and save variables*/
   anglim () ;
   shkang = thetmax ;
   compr () ;
   adet  = ang     ;
   pdet  = prat    ;
   tdet  = trat    ;
   rdet  = rhorat  ;
   ptdet = ptrat   ;
   m2det = mach2   ;
/*      find the limiting value of the other variables and save them */
   shkang = 90.0   ;
   compr () ;
   pmax  = prat    ;
   tmax  = trat    ;
   rmax  = rhorat  ;
   ptmin = ptrat   ;
   m2min = mach2   ;

    switch (inopt2) {
        case 1: {                             /* ramp angle  given */
             ang = fl_getr(ino2) ;
             if (ang >= adet) {erflo=1; fl_set_input(outo,"Detached"); return ;}
             if (ang <= 0.0) {erflo=1; fl_set_input(outo,"Expansion"); return ;}
             get_shang_ang (ang*convdr) ;
             break ;
         }
        case 2: {                             /* shock angle given */
             shkang = fl_getr(ino2) ;
             if (shkang >= thetmax) {
                 erflo=1; fl_set_input(outo,"Detached"); return ;}
             if (shkang <= asin(1.0/mach1)/convdr) {
                 erflo=1; fl_set_input(outo,"Shock < M ang"); return ;}
             break ;
         }
        case 3: {                             /* total press ratio given */
             ptrat = fl_getr(ino2) ;
             if (ptrat <= ptmin) {erflo=1; fl_set_input(outo,"Detached");return;}
             if (ptrat > 1.0) {erflo=1; fl_set_input(outo,"pt2/pt1>1"); return; }
             msave = mach1 ;
             get_mach_ptrat() ;
             shkang = asin(mach1/msave) / convdr ;
             mach1 = msave ;
             break ;
         }
        case 4: {                             /* static pressure ratio given */
             prat = fl_getr(ino2) ;
             if (prat >= pmax){erflo=1; fl_set_input(outo,"Detached"); return; }
             if (prat <= 1.0){erflo=1; fl_set_input(outo,"Expansion"); return; }
             msint = sqrt((gp1 * prat + gm1)/(2.0*gama)) ;
             shkang = asin(msint/mach1) / convdr ;
             break ;
         }
        case 5: {                             /* temperature ratio given */
             trat = fl_getr(ino2) ;
             if (trat >= tmax){erflo=1; fl_set_input(outo,"Detached"); return; }
             if (trat <= 1.0){erflo=1; fl_set_input(outo,"Expansion"); return; }
             msave = mach1 ;
             get_mach_trat() ;
             shkang = asin(mach1/msave) / convdr ;
             mach1 = msave ;
             break ;
         }
        case 6: {                             /* density ratio given */
             rhorat = fl_getr(ino2) ;
             if (rhorat >=rmax){erflo=1; fl_set_input(outo,"Detached"); return; }
             if (rhorat <= 1.0){erflo=1;fl_set_input(outo,"Expansion"); return; }
             msint = sqrt(2.0*rhorat/(gp1 - gm1*rhorat)) ;
             shkang = asin(msint/mach1) / convdr ;
             break ;
         }
        case 7: {                          /* downstream Mach number given */
             mach2 = fl_getr(ino2) ;
             if (mach2 >= mach1){erflo=1;fl_set_input(outo,"Expansion");return; }
             if (mach2 <= m2min){erflo=1;fl_set_input(outo,"Detached");return; }
             get_shang_mach2 () ;
             break ;
         }
     }
     compr () ;               
     return ;
}

void get_shang_ang (float delr)            /* oblique shock problem      */
{                                          /*  get shock angle if Mach 1 */
   float sints,cotd,mach1s ;               /*  and ramp angle given      */
   float theto,thetn,delo,deln,deriv ;

   mach1s = mach1*mach1 ;

   theto = asin(1.0/mach1) ;             /* iterate for shock angle */
   delo = 0.0 ;
   thetn = theto + 3.0 * convdr ;
   while (fabs(delr - delo) > .0001) {
         sints = pow(sin(thetn),2.0) ;
         cotd = tan(thetn)*((gp1*mach1s)/         /* EQ. 138 */
                (2.0*(mach1s * sints - 1.0))-1.0);
         deln = atan(1.0/cotd) ;
         deriv = (deln-delo)/(thetn-theto) ;
         delo = deln ;
         theto = thetn ;
         thetn = theto + (delr-delo)/deriv ;
    }

   shkang = theto / convdr ;
   return;
}

void get_shang_mach2 ()                    /* oblique shock problem      */
{                                          /*  get shock angle if Mach 1 */
   float mach1s,mach1f;                    /*  and Mach 2 given          */
   float sintso,sintsn,mach2o,mach2n,deriv ;

   mach1s = mach1*mach1 ;
   mach1f = mach1s*mach1s ;

   sintso = 1.0/mach1s ;
   mach2o = mach1 ;
   sintsn = sintso + .05 ;
   while (fabs(mach2 - mach2o) > .0001) {
         mach2n = sqrt(((mach1f * sintsn * pow(gp1,2.0)) +      /* EQ. 132 */
            (-4.0 * (mach1s*sintsn-1.0)*(gama*mach1s*sintsn+1.0))) /
            ((2.0*gama*mach1s*sintsn - gm1)*(gm1*mach1s*sintsn + 2.0))) ;
         deriv = (mach2n-mach2o)/(sintsn-sintso) ;
         mach2o = mach2n ;
         sintso = sintsn ;
         sintsn = sintso + (mach2-mach2o)/deriv ;
    }

   shkang = asin(sqrt(sintso)) / convdr ;
   return;
}

void anglim()                     /* determine the max attach shock angle  */
{                                        /*  for a given  Mach number      */
   float a1,b1,c1,sints,machs,machf ;
   float cotd,delmr ;

   machs = mach1*mach1 ;
   machf = machs*machs ;

   a1 = 2.0 * gama * machf ;         /* max from my derivation for large ang*/
   b1 = 4.0 * machs - gp1 * machf ;
   c1 = -(2.0 + gp1 * machs) ;
   sints = (-b1 + sqrt(pow(b1,2.0)-4.0*a1*c1))/(2.0*a1) ;
   thetmax = asin(sqrt(sints)) /convdr ;

   cotd = tan(convdr*thetmax)*(((gama+1.0)*machs)/    /* EQ. 138 */
         (2.0*(machs * sints - 1.0))-1.0);
   delmr = atan(1.0/cotd) ;
   delmax = delmr /convdr ;

   return ;
}

void compr ()                               /*  oblique shock relations  */
{                                           /*  equations from NACA 1135 */
   float mach2s,sints,cotd ;                /* Mach 1 and shock angle given */
   float mach1s,mach1f ;

   mach1s = mach1*mach1 ;
   mach1f = mach1s*mach1s ;
   sints = pow(sin(shkang*convdr),2.0) ;

   prat = (2.0*gama*mach1s*sints - gm1)/gp1 ;                  /*  EQ. 128 */
   rhorat = (gp1*mach1s*sints)/(gm1*mach1s*sints + 2.0) ;      /*  EQ. 129 */
   trat = (2.0*gama*mach1s*sints - gm1) * (gm1*mach1s*sints + 2.0)
                 /(mach1s*sints*pow(gp1,2.0)) ;                /* EQ. 130 */
   mach2s = ((mach1f * sints * pow(gp1,2.0)) +
            (-4.0 * (mach1s*sints-1.0)*(gama*mach1s*sints+1.0))) /
     ((2.0*gama*mach1s*sints - gm1)*(gm1*mach1s*sints + 2.0)) ; /* EQ. 132 */
   mach2 = sqrt(mach2s) ;
   cotd = tan(shkang*convdr)*((gp1*mach1s)/                     /* EQ. 138 */
            (2.0*(mach1s * sints - 1.0))-1.0);
   ang = (atan(1.0/cotd))/convdr ;
   ptrat = (pow(((gp1*mach1s*sints)/(gm1*mach1s*sints + 2.0)),(gama/gm1)))
       * pow((gp1/(2.0*gama*mach1s*sints - gm1)),(1.0/gm1)) ;   /* EQ. 142 */

   return;
}

/*      ------------    Centered Expansion Calculator  ------- */

void exp_swtch (FL_OBJECT *obj, long arg)
/*   switching routine for Centered Expansion - uses radio button*/
{
  void anlye () ;

  gama = fl_getr(game) ;          /* get gamma */
  anlye () ;
  if (erflo != 1) {
     switch (arg)
     {                /* output */
      case 1: fl_putr (oute,mach1,4) ; break ;  /* upstream Mach */
      case 2: fl_putr (oute,mach2,4) ; break ;  /* downstream Mach */
      case 3: fl_putr (oute,ang,2) ;   break ;  /* ramp angle */
      case 4: fl_putr (oute,shkang,2); break ;  /* upstream Mach angle */
      case 5: fl_putr (oute,expang,2); break ;  /* downstream Mach angle */
      case 6: fl_putr (oute,prat,4) ;  break ;  /* static pressure */
      case 7: fl_putr (oute,trat,4) ;  break ;  /* temperature */
      case 8: fl_putr (oute,rhorat,4); break ;  /* density */
     }
  }
}

void create_eform()
{
  FL_OBJECT *obj;

  eform = fl_bgn_form(FL_UP_BOX,290.0,220.0);

               /* Input */
  obj = fl_add_box(FL_NO_BOX,5.0,190.0,90.0,25.0,"Mach 1");
  ine1 = fl_add_input(FL_NORMAL_INPUT,5.,160.,90.,25.,"") ;
     fl_set_object_color(ine1,7,7);
     fl_set_object_lsize(ine1,FL_MEDIUM_SIZE); 
  che2 = fl_add_choice(FL_NORMAL_CHOICE,100.0,190.0,85.0,25.0,"");
     fl_set_object_boxtype(che2,FL_RSHADOW_BOX);
     fl_addto_choice(che2,"Ramp");
     fl_addto_choice(che2,"Mach 2");
     fl_addto_choice(che2,"p2/p1");
     fl_addto_choice(che2,"T2/T1");
     fl_addto_choice(che2,"r2/r1");
     fl_addto_choice(che2,"M2-ang");
  ine2 = fl_add_input(FL_NORMAL_INPUT,100.,160.,85.,25.,"") ;
     fl_set_object_color(ine2,7,7);
     fl_set_object_lsize(ine2,FL_MEDIUM_SIZE); 
             /* gamma window */
    obj = fl_add_box(FL_NO_BOX,195.0,190.0,90.0,25.0,"g");
       fl_set_object_lstyle(obj,15) ;
       fl_set_object_lsize(obj,FL_NORMAL_SIZE);
    game = obj = fl_add_input(FL_NORMAL_INPUT,195.,160.,90.,25.,"") ;
       fl_set_object_color(obj,7,7);
       fl_set_object_lsize(obj,FL_MEDIUM_SIZE);

             /* Choice of Outputs */
  obj = fl_add_box(FL_FLAT_BOX,5.0,30.0,280.0,125.0,"");
    fl_set_object_color(obj,7,7);
  obj = fl_bgn_group ();
     obj = fl_add_lightbutton(FL_RADIO_BUTTON,5.0,125.0,90.0,25.0,"Mach 1");
       fl_set_object_callback(obj,exp_swtch,1);
     obj = fl_add_lightbutton(FL_RADIO_BUTTON,100.0,125.0,90.0,25.0,"Mach 2");
       fl_set_object_callback(obj,exp_swtch,2);
     obj = fl_add_lightbutton(FL_RADIO_BUTTON,195.0,125.0,90.0,25.0,"Ramp");
       fl_set_object_callback(obj,exp_swtch,3);
     obj = fl_add_lightbutton(FL_RADIO_BUTTON,5.0,95.0,90.0,25.0,"M1-ang");
       fl_set_object_callback(obj,exp_swtch,4);
     obj = fl_add_lightbutton(FL_RADIO_BUTTON,100.0,95.0,90.0,25.0,"M2-ang");
       fl_set_object_callback(obj,exp_swtch,5);
     obj = fl_add_lightbutton(FL_RADIO_BUTTON,195.0,95.0,90.0,25.0,"p2/p1");
       fl_set_object_callback(obj,exp_swtch,6);
     obj = fl_add_lightbutton(FL_RADIO_BUTTON,55.0,65.0,90.0,25.0,"T2/T1");
       fl_set_object_callback(obj,exp_swtch,7);
     obj = fl_add_lightbutton(FL_RADIO_BUTTON,150.0,65.0,90.0,25.0,"r2/r1");
       fl_set_object_callback(obj,exp_swtch,8);
             /* Output window */
  oute = fl_add_input(FL_NORMAL_INPUT,100.,35.,100.,25.,"") ;
       fl_set_input_color(oute,0,7);
       fl_set_object_lsize(oute,FL_MEDIUM_SIZE);

              /* Exit Button */
 obj = fl_add_button(FL_NORMAL_BUTTON,100.0,5.0,90.0,25.0,"Exit");
         fl_set_object_callback(obj,clos_swtch,2); 
         fl_set_object_color(obj,1,1);
         fl_set_object_lcol(obj,7);
  fl_end_form();
}

void anlye ()                           /* analysis for centered expansion */
{
    float m1,p1,t1,r1,nu1;
    void isen () ;
    void get_machpm ();

    erflo = 0 ;
    if (gama <1.01) gama = 1.01 ;
    gm1 = gama - 1.0 ;
    gp1 = gama + 1.0;
    inopt2 = fl_get_choice(che2);
    mach1 = fl_getr(ine1);          /* get mach number*/
    if (mach1 < 1.0) {erflo=1; fl_set_input(oute,"M < 1"); return; }

    isen ();
    m1 = mach1 ;
    shkang = mu;
    p1 = prat;
    t1 = trat;
    r1 = rhorat ;
    nu1 = nu ;

    switch (inopt2) {
        case 1: {                             /* wedge angle given */
             ang = fl_getr(ine2) ;
/*   if (ang <= 0.0) {erflo=1; fl_set_input(oute,"neg angle"); return; } */
             nu = ang + nu1;
             get_machpm ();
             break;
         }
         case 2: {                            /* downstream Mach given */
             mach2 = fl_getr(ine2);
             if (mach2 <= mach1) {erflo=1; fl_set_input(oute,"M2 < M1"); return;}
             mach1 = mach2 ;
             break;
         }
         case 3: {                               /* static pressure given */
             p2p1 = fl_getr(ine2);
             if (p2p1 >= 1.0) {erflo=1; fl_set_input(oute,"p2 > p1"); return; }
             mach1 = sqrt(2.0*(pow(p2p1*p1,-gm1/gama)-1.0)/gm1) ;
             break ;
         }
        case 4: {                             /* temperature ratio given */
             t2t1 = fl_getr(ine2);
             if (t2t1 >= 1.0) {erflo=1; fl_set_input(oute,"T2 > T1"); return; }
             mach1 = sqrt(2.0*(pow(t2t1*t1,-1.0) - 1.0)/gm1) ;
             break ;
         }
        case 5: {                             /* density ratio given */
             r2r1 = fl_getr(ine2);
             if (r2r1 >= 1.0) {erflo=1; fl_set_input(oute,"r2 > r1"); return; }
             mach1 = sqrt(2.0*(pow(r2r1*r1,-gm1) - 1.0)/gm1) ;
             break ;
         }
        case 6: {                             /* downstream mach angle given */
             ang = fl_getr(ine2) ;
             mach1 = 1.0 / sin(ang*convdr) ;
             break ;
         }
     }
     isen ();
     ang = nu - nu1 ;
 
     mach2 = mach1 ;
     mach1 = m1 ;
     expang = mu ;
     prat = prat / p1 ;
     trat = trat / t1 ;
     rhorat = rhorat / r1 ;

     return ;
}

/*      ------------    Rayleigh Flow Calculator  ------- */

void ray_swtch (FL_OBJECT *obj, long arg)
/*   callback routine for Rayleigh Flow */
{
  void anlyr () ;

  gama = fl_getr(gamr) ;          /* get gamma */
  anlyr () ;
  if (erflo != 1) {
     switch (arg)
     {                /* output */
      case 1: fl_putr (outr,mach1,3) ; break ;  /*  Mach */
      case 2: fl_putr (outr,ttrat,3); break ;    /* temperature ratio */
      case 3: fl_putr (outr,ptrat,3) ; break ;   /* pressure ratio */
      case 4: fl_putr (outr,mach2,3) ;  break ; /* downstream mach */
      case 5: fl_putr (outr,prat,3) ;  break ; /* pressure ratio */
      case 6: fl_putr (outr,trat,3) ;  break ; /* temp ratio */
      case 7: fl_putr (outr,vrat,3) ;  break ; /* temp ratio */
      case 8: fl_putr (outr,rhorat,3) ;  break ; /* temp ratio */
     }
  }
}

void create_rform()
{
  FL_OBJECT *obj;

  rform = fl_bgn_form(FL_UP_BOX,290.0,220.0);

               /* Input */
  chr1 = fl_add_choice(FL_NORMAL_CHOICE,5.0,190.0,90.0,25.0,"");
     fl_set_object_boxtype(chr1,FL_RSHADOW_BOX);
     fl_addto_choice(chr1,"Mach 1");
     fl_addto_choice(chr1,"Mach 2");
  inr1 = fl_add_input(FL_NORMAL_INPUT,5.,160.,90.,25.,"") ;
     fl_set_object_color(inr1,7,7);
     fl_set_object_lsize(inr1,FL_MEDIUM_SIZE);
  chr2 = fl_add_choice(FL_NORMAL_CHOICE,100.0,190.0,85.0,25.0,"");
     fl_set_object_boxtype(chr2,FL_RSHADOW_BOX);
     fl_addto_choice(chr2,"Mach 2");
     fl_addto_choice(chr2,"T2/T1");
     fl_addto_choice(chr2,"p2/p1");
     fl_addto_choice(chr2,"Tt2/Tt1");
     fl_addto_choice(chr2,"pt2/pt1");
     fl_addto_choice(chr2,"V2/V1");
     fl_addto_choice(chr2,"rho2/rho1");
  inr2 = fl_add_input(FL_NORMAL_INPUT,100.,160.,85.,25.,"") ;
     fl_set_object_color(inr2,7,7);
     fl_set_object_lsize(inr2,FL_MEDIUM_SIZE);
             /* gamma window */
  obj = fl_add_box(FL_NO_BOX,195.0,190.0,90.0,25.0,"g");
       fl_set_object_lstyle(obj,15) ;
       fl_set_object_lsize(obj,FL_NORMAL_SIZE);
  gamr = obj = fl_add_input(FL_NORMAL_INPUT,195.,160.,90.,25.,"") ;
     fl_set_object_color(obj,7,7);
     fl_set_object_lsize(obj,FL_MEDIUM_SIZE);

             /* Choice of Outputs */
  obj = fl_add_box(FL_FLAT_BOX,5.0,30.0,280.0,125.0,"");
           fl_set_object_color(obj,7,7);
     obj = fl_add_lightbutton(FL_RADIO_BUTTON,5.0,125.0,90.0,25.0,"Mach 1");
       fl_set_object_callback(obj,ray_swtch,1);
     obj = fl_add_lightbutton(FL_RADIO_BUTTON,100.0,125.0,90.0,25.0,"Mach 2");
       fl_set_object_callback(obj,ray_swtch,4);
     obj = fl_add_lightbutton(FL_RADIO_BUTTON,195.0,125.0,90.0,25.0,"pt2/pt1");
       fl_set_object_callback(obj,ray_swtch,3);
     obj = fl_add_lightbutton(FL_RADIO_BUTTON,5.0,95.0,90.0,25.0,"Tt2/Tt1");
       fl_set_object_callback(obj,ray_swtch,2);
     obj = fl_add_lightbutton(FL_RADIO_BUTTON,100.0,95.0,90.0,25.0,"p2/p1");
       fl_set_object_callback(obj,ray_swtch,5);
     obj = fl_add_lightbutton(FL_RADIO_BUTTON,195.0,95.0,90.0,25.0,"T2/T1" );
       fl_set_object_callback(obj,ray_swtch,6);
     obj = fl_add_lightbutton(FL_RADIO_BUTTON,5.0,65.0,90.0,25.0,"V2/V1" );
       fl_set_object_callback(obj,ray_swtch,7);
     obj = fl_add_lightbutton(FL_RADIO_BUTTON,100.0,65.0,90.0,25.0,"rho2/rho1");
       fl_set_object_callback(obj,ray_swtch,8);
             /* Output window */
  outr = fl_add_input(FL_NORMAL_INPUT,100.,35.,100.,25.,"") ;
       fl_set_input_color(outr,0,7);
       fl_set_object_lsize(outr,FL_MEDIUM_SIZE);

              /* Exit Button */
  obj = fl_add_button(FL_NORMAL_BUTTON,100.0,5.0,90.0,25.0,"Exit");
         fl_set_object_callback(obj,clos_swtch,7); 
         fl_set_object_color(obj,1,1);
         fl_set_object_lcol(obj,7);
  fl_end_form();
}

void anlyr ()                           /* analysis for rayleigh flow */
{                                       /*  reworked - Jan 30 - 96    */
   float mread,ratio,rstar,machs,gp1,gm1 ;
   float mold,mnew,rold,rnew,deriv ;
   float get_pray(float,float);
   float get_tray(float,float);
   float get_ptray(float,float);
   float get_ttray(float,float);
   float get_vray(float,float);
   void rayleigh(float,float,float) ;

   erflo = 0 ;
   gp1 = gama + 1.0 ;
   gm1 = gama - 1.0 ;
   inopt1 = fl_get_choice(chr1);
   inopt2 = fl_get_choice(chr2);

   mread = fl_getr(inr1) ;
   ratio = fl_getr(inr2) ;

   if (inopt1 == 1) mach1 = mread ;
   else {
          mach2 = mread ;
          if(inopt2 != 1) ratio = 1.0/ratio ;
   }

   switch (inopt2) {
       case 1: {                /* mach2 given */
               if (inopt1 == 2) {   /* mach2 given */
                   erflo =1 ;
                   fl_set_input(outr,"ERROR");
                   return ;
               }
               else {              /* mach1 given */
                   mach2 = ratio ;
               }
               if (mach1 < 1.0 && mach2 > 1.0) {
                   erflo =1 ;
                   fl_set_input(outr,"M2 > 1");
                   return ;
               }
               if (mach1 > 1.0 && mach2 < 1.0) {
                   erflo =1 ;
                   fl_set_input(outr,"M2 < 1");
                   return ;
               }
               break;
       }
       case 2: {                /* T2/T1 given */
               rstar = get_tray(mread,gama)*ratio ;
               if (rstar > get_tray(1.0/sqrt(gama),gama)) {
                   erflo =1 ;
                   fl_set_input(outr,"Choked");
                   return ;
               }
               rold = 1.0 ;
               mold = 1.0 ;
               mnew = mread ;
               while (fabs(rstar - rold) > .0001) {
                     rnew = get_tray(mnew,gama) ;
                     deriv = (rnew-rold)/(mnew-mold) ;
                     rold = rnew ;
                     mold = mnew ;
                     mnew = mold + (rstar - rold)/deriv ;
               }
               if (inopt1 == 2) mach1 = mold ;
               else mach2 = mold ;
               break;
       }
       case 3: {                /* p2/p1 given */
               rstar = get_pray(mread,gama)*ratio ;
               if (mread < 1.0 && rstar > gama+1.0 ) {
                   erflo =1 ;
                   fl_set_input(outr,"Stagnate");
                   return ;
               }
               if (mread > 1.0 && rstar > 1.0) {
                      erflo =1 ;
                      fl_set_input(outr,"Choked");
                      return ;
               }
               if (mread < 1.0 && rstar < 1.0) {
                      erflo =1 ;
                      fl_set_input(outr,"Choked");
                      return ;
               }
               machs = ((gp1/rstar)-1.0)/gama ;
               if(inopt1 == 2) mach1 = sqrt(machs) ;
               else mach2 = sqrt(machs) ;
               break;
       }
       case 4: {                /* Tt2/Tt1 given */
               rstar = get_ttray(mread,gama)*ratio ;
               if (rstar > 1.0 ) {
                   erflo =1 ;
                   fl_set_input(outr,"Choked");
                   return ;
               }
               if (mread > 1.0 && rstar < 0.5) {
                      erflo =1 ;
                      fl_set_input(outr,"Impossible");
                      return ;
               }
               rold = 1.0 ;
               mold = 1.0 ;
               mnew = mread ;
               while (fabs(rstar - rold) > .0001) {
                     rnew = get_ttray(mnew,gama) ;
                     deriv = (rnew-rold)/(mnew-mold) ;
                     rold = rnew ;
                     mold = mnew ;
                     mnew = mold + (rstar - rold)/deriv ;
               }
               if (inopt1 == 2) mach1 = mold ;
               else mach2 = mold ;
               break;
       }
       case 5: {                /* pt2/pt1 given */
               rstar = get_ptray(mread,gama)*ratio ;
               if (rstar < 1.0 ) {
                   erflo =1 ;
                   fl_set_input(outr,"Choked");
                   return ;
               }
               if (mread < 1.0 && rstar > get_ptray(0.0,gama) ) {
                   erflo =1 ;
                   fl_set_input(outr,"Stagnate");
                   return ;
               }
               rold = 1.0 ;
               mold = 1.0 ;
               mnew = mread ;
               while (fabs(rstar - rold) > .0001) {
                     rnew = get_ptray(mnew,gama) ;
                     deriv = (rnew-rold)/(mnew-mold) ;
                     rold = rnew ;
                     mold = mnew ;
                     mnew = mold + (rstar - rold)/deriv ;
               }
               if (inopt1 == 2) mach1 = mold ;
               else mach2 = mold ;
               break;
       }
       case 6: {                /* V2/V1 given */
               rstar = get_vray(mread,gama)*ratio ;
               if (mread > 1.0 && rstar < 1.0) {
                      erflo =1 ;
                      fl_set_input(outr,"Choked");
                      return ;
               }
               if (mread < 1.0 && rstar > 1.0) {
                      erflo =1 ;
                      fl_set_input(outr,"Choked");
                      return ;
               }
               if (mread > 1.0 && rstar > 1.7) {
                      erflo =1 ;
                      fl_set_input(outr,"Impossible");
                      return ;
               }
               machs = rstar/(gp1-gama*rstar) ;
               if (inopt1 == 2)  mach1 = sqrt(machs) ;
               else mach2 = sqrt(machs) ;
               break;
       }
       case 7: {                /* rho2/rho1 given */
               ratio = 1.0/ratio ;
               rstar = get_vray(mread,gama)*ratio ;
               if (mread > 1.0 && rstar < 1.0) {
                       erflo =1 ;
                       fl_set_input(outr,"Choked");
                       return ;
               }
               if (mread < 1.0 && rstar > 1.0) {
                       erflo =1 ;
                       fl_set_input(outr,"Choked");
                       return ;
               }
               if (mread > 1.0 && rstar > 1.7) {
                      erflo =1 ;
                      fl_set_input(outr,"Impossible");
                      return ;
               }
               machs = rstar/(gp1-gama*rstar) ;
               if (inopt1 == 2)  mach1 = sqrt(machs) ;
               else mach2 = sqrt(machs) ;
               break;
       }
   }

   rayleigh(mach1,mach2,gama) ;  
   return ;
}

float get_pray(float machin, float gama)
/* Utility to get the rayleigh pressure ratio given the mach number */
{
    float number;
    number = (gama+1.0)/(1.0+gama*machin*machin) ;

    return(number) ;
}

float get_ptray(float machin, float gama)
/* Utility to get the rayleigh total pressure ratio given the mach number */
{
    float number,fac1,gm1,gp1,mach1s;
    mach1s = machin*machin ;
    gp1 = gama + 1.0 ;
    gm1 = gama - 1.0 ;
    fac1 = 2.0*(1.0+.5*gm1*mach1s)/gp1 ;
    number = pow(fac1,(gama/gm1))*gp1/(1.0+gama*machin*machin) ;

    return(number) ;
}

float get_tray(float machin, float gama)
/* Utility to get the rayleigh temperature ratio given the mach number */
{
    float number,mach1s,gp1;
    mach1s = machin*machin ;
    gp1 = gama + 1.0 ;
    number = gp1*gp1*mach1s/pow((1.0+gama*machin*machin),2.0) ;

    return(number) ;
}

float get_ttray(float machin, float gama)
/* Utility to get the rayleigh total temp ratio given the mach number */
{
    float number,fac1,gm1,gp1,mach1s;
    mach1s = machin*machin ;
    gp1 = gama + 1.0 ;
    gm1 = gama - 1.0 ;
    fac1 = 2.0*(1.0+.5*gm1*mach1s)*gp1*mach1s ;
    number = fac1/pow((1.0+gama*machin*machin),2.0) ;

    return(number) ;
}

float get_vray(float machin, float gama)
/* Utility to get the rayleigh velocity ratio given the mach number */
{
    float number;
    number = (gama+1.0)*machin*machin/(1.0+gama*machin*machin) ;

    return(number) ;
}

void rayleigh (float machin, float machout, float gama)
{                           /* get rayleigh flow ratios */
   float get_pray(float,float);
   float get_tray(float,float);
   float get_ptray(float,float);
   float get_ttray(float,float);
   float get_vray(float,float);
 
   prat  = get_pray(machout,gama) /get_pray(machin,gama) ;
   ptrat = get_ptray(machout,gama)/get_ptray(machin,gama) ;
   trat  = get_tray(machout,gama) /get_tray(machin,gama) ;
   ttrat = get_ttray(machout,gama)/get_ttray(machin,gama) ;
   vrat  = get_vray(machout,gama) /get_vray(machin,gama) ;
   rhorat = 1.0/vrat ;
}
   
/*      ------------  Fanno  Flow Calculator  ------- */

void fan_swtch (FL_OBJECT *obj, long arg)
/*   callback routine for Fanno Flow */
{
  void anlyf () ;

  gama = fl_getr(gamf) ;          /* get gamma */
  anlyf () ;
  if (erflo != 1) {
     switch (arg)
     {                /* output */
      case 1: fl_putr (outf,mach1,3) ; break ;  /*  Mach */
      case 2: fl_putr (outf,frat,3)  ; break ;  /* impule ratio */
      case 3: fl_putr (outf,ptrat,3) ; break ;   /* pressure ratio */
      case 4: fl_putr (outf,mach2,3) ;  break ;  /* downstream mach */
      case 5: fl_putr (outf,prat,3) ;  break ; /* pressure ratio */
      case 6: fl_putr (outf,trat,3) ;  break ; /* temp ratio */
      case 7: fl_putr (outf,vrat,3) ;  break ; /* temp ratio */
      case 8: fl_putr (outf,rhorat,3) ;  break ; /* temp ratio */
      case 9: fl_putr (outf,leng,3) ;  break ; /* length of duct */
     }
  }
}

void create_fform()
{
  FL_OBJECT *obj;

  fform = fl_bgn_form(FL_UP_BOX,290.0,220.0);

               /* Input */
  chf1 = fl_add_choice(FL_NORMAL_CHOICE,5.0,190.0,90.0,25.0,"");
     fl_set_object_boxtype(chf1,FL_RSHADOW_BOX);
     fl_addto_choice(chf1,"Mach 1");
     fl_addto_choice(chf1,"Mach 2");
  inf1 = fl_add_input(FL_NORMAL_INPUT,5.,160.,90.,25.,"") ;
     fl_set_object_color(inf1,7,7);
     fl_set_object_lsize(inf1,FL_MEDIUM_SIZE);
  chf2 = fl_add_choice(FL_NORMAL_CHOICE,100.0,190.0,85.0,25.0,"");
     fl_set_object_boxtype(chf2,FL_RSHADOW_BOX);
     fl_addto_choice(chf2,"Mach 2");
     fl_addto_choice(chf2,"T2/T1");
     fl_addto_choice(chf2,"p2/p1");
     fl_addto_choice(chf2,"F2/F1");
     fl_addto_choice(chf2,"pt2/pt1");
     fl_addto_choice(chf2,"V2/V1");
     fl_addto_choice(chf2,"rho2/rho1");
     fl_addto_choice(chf2,"4fL/D");
  inf2 = fl_add_input(FL_NORMAL_INPUT,100.,160.,85.,25.,"") ;
     fl_set_object_color(inf2,7,7);
     fl_set_object_lsize(inf2,FL_MEDIUM_SIZE);
             /* gamma window */
  obj = fl_add_box(FL_NO_BOX,195.0,190.0,90.0,25.0,"g");
       fl_set_object_lstyle(obj,15) ;
       fl_set_object_lsize(obj,FL_NORMAL_SIZE);
  gamf = obj = fl_add_input(FL_NORMAL_INPUT,195.,160.,90.,25.,"") ;
     fl_set_object_color(obj,7,7);
     fl_set_object_lsize(obj,FL_MEDIUM_SIZE);

             /* Choice of Outputs */
  obj = fl_add_box(FL_FLAT_BOX,5.0,30.0,280.0,125.0,"");
           fl_set_object_color(obj,7,7);
     obj = fl_add_lightbutton(FL_RADIO_BUTTON,5.0,125.0,90.0,25.0,"Mach 1");
       fl_set_object_callback(obj,fan_swtch,1);
     obj = fl_add_lightbutton(FL_RADIO_BUTTON,100.0,125.0,90.0,25.0,"Mach 2");
       fl_set_object_callback(obj,fan_swtch,4);
     obj = fl_add_lightbutton(FL_RADIO_BUTTON,195.0,125.0,90.0,25.0,"pt2/pt1");
       fl_set_object_callback(obj,fan_swtch,3);
     obj = fl_add_lightbutton(FL_RADIO_BUTTON,5.0,95.0,90.0,25.0,"F2/F1");
       fl_set_object_callback(obj,fan_swtch,2);
     obj = fl_add_lightbutton(FL_RADIO_BUTTON,100.0,95.0,90.0,25.0,"p2/p1");
       fl_set_object_callback(obj,fan_swtch,5);
     obj = fl_add_lightbutton(FL_RADIO_BUTTON,195.0,95.0,90.0,25.0,"T2/T1" );
       fl_set_object_callback(obj,fan_swtch,6);
     obj = fl_add_lightbutton(FL_RADIO_BUTTON,5.0,65.0,90.0,25.0,"V2/V1" );
       fl_set_object_callback(obj,fan_swtch,7);
     obj = fl_add_lightbutton(FL_RADIO_BUTTON,100.0,65.0,90.0,25.0,"rho2/rho1");
       fl_set_object_callback(obj,fan_swtch,8);
     obj = fl_add_lightbutton(FL_RADIO_BUTTON,195.0,65.0,90.0,25.0,"4fL/D");
       fl_set_object_callback(obj,fan_swtch,9);
             /* Output window */
  outf = fl_add_input(FL_NORMAL_INPUT,100.,35.,100.,25.,"") ;
       fl_set_input_color(outf,0,7);
       fl_set_object_lsize(outf,FL_MEDIUM_SIZE);

              /* Exit Button */
  obj = fl_add_button(FL_NORMAL_BUTTON,100.0,5.0,90.0,25.0,"Exit");
         fl_set_object_callback(obj,clos_swtch,9); 
         fl_set_object_color(obj,1,1);
         fl_set_object_lcol(obj,7);

  fl_end_form();
}

void anlyf ()                           /* analysis for fanno flow */
{                                  
   float mread,ratio,rstar,machs,gp1,gm1 ;
   float mold,mnew,rold,rnew,deriv ;
   int iter ;
   float get_pfan(float,float);
   float get_vfan(float,float);
   float get_tfan(float,float);
   float get_ptfan(float,float);
   float get_ffan(float,float);
   float get_lfan(float,float);
   void fanno(float,float,float) ;

   erflo = 0 ;
   gp1 = gama + 1.0 ;
   gm1 = gama - 1.0 ;
   inopt1 = fl_get_choice(chf1);
   inopt2 = fl_get_choice(chf2);

   mread = fl_getr(inf1) ;
   ratio = fl_getr(inf2) ;

   if (inopt1 == 1) mach1 = mread ;
   else mach2 = mread ;

   switch (inopt2) {
       case 1: {                /* mach2 given */
               if (inopt1 == 2) {   /* mach2 given */
                   erflo =1 ;
                   fl_set_input(outf,"ERROR");
                   return ;
               }
               else {              /* mach1 given */
                   mach2 = ratio ;
               }
               if (mach1 < 1.0 && mach2 < mach1) {
                   erflo =1 ;
                   fl_set_input(outf,"M2<M1");
                   return ;
               }
               if (mach1 > 1.0 && mach2 > mach1) {
                   erflo =1 ;
                   fl_set_input(outf,"M2>M1");
                   return ;
               }
               if (mach1 < 1.0 && mach2 > 1.0) {
                   erflo =1 ;
                   fl_set_input(outf,"M2 > 1");
                   return ;
               }
               if (mach1 > 1.0 && mach2 < 1.0) {
                   erflo =1 ;
                   fl_set_input(outf,"M2 < 1");
                   return ;
               }
               break;
       }
       case 2: {                /* T2/T1 given */
               if (mread < 1.0 && ratio > 1.0) {
                   erflo =1 ;
                   fl_set_input(outf,"T2>T1");
                   return ;
               }
               if (mread > 1.0 && ratio < 1.0) {
                   erflo =1 ;
                   fl_set_input(outf,"T2<T1");
                   return ;
               }
               if(inopt1 == 2) ratio = 1.0/ratio ;
               rstar = get_tfan(mread,gama)*ratio ;
               if (mread < 1.0 && rstar > gp1/2.0) {
                   erflo =1 ;
                   fl_set_input(outf,"Stagnate");
                   return ;
               }
               if (mread < 1.0 && rstar < 1.0) {
                   erflo =1 ;
                   fl_set_input(outf,"Choked");
                   return ;
               }
               if (mread > 1.0 && rstar > 1.0) {
                   erflo =1 ;
                   fl_set_input(outf,"Choked");
                   return ;
               }
               machs = (gp1-2.0*rstar)/(gm1*rstar);

               if(inopt1 == 2) mach1 = sqrt(machs) ;
               else mach2 = sqrt(machs) ;
               break;
       }
       case 3: {                /* p2/p1 given */
               if (mread < 1.0 && ratio > 1.0) {
                   erflo =1 ;
                   fl_set_input(outf,"p2>p1");
                   return ;
               }
               if (mread > 1.0 && ratio < 1.0) {
                   erflo =1 ;
                   fl_set_input(outf,"p2<p1");
                   return ;
               }
               if(inopt1 == 2) ratio = 1.0/ratio ;
               rstar = get_pfan(mread,gama)*ratio ;
               if (mread < 1.0 && rstar < 1.0) {
                   erflo =1 ;
                   fl_set_input(outf,"Choked");
                   return ;
               }
               if (mread > 1.0 && rstar > 1.0) {
                   erflo =1 ;
                   fl_set_input(outf,"Choked");
                   return ;
               }
               machs = (-1.0+sqrt(1.0+gp1*gm1/(rstar*rstar)))/gm1 ;

               if(inopt1 == 2) mach1 = sqrt(machs) ;
               else mach2 = sqrt(machs) ;
               break;
       }
       case 4: {                /* F2/F1 given */
               if (ratio > 1.0) {
                   erflo =1 ;
                   fl_set_input(outf,"F2>F1");
                   return ;
               }
               if(inopt1 == 2) ratio = 1.0/ratio ;
               rstar = get_ffan(mread,gama)*ratio ;
               if (rstar <= 1.0) {
                   erflo =1 ;
                   fl_set_input(outf,"Choked");
                   return ;
               }
               if (mread > 1.0 && rstar > gama) {
                   erflo =1 ;
                   fl_set_input(outf,"F/F*>gam");
                   return ;
               }
               rold = 1.0 ;
               mold = 1.0 ;
               mnew = mread ;
               while (fabs(rstar - rold) > .0001) {
                     rnew = get_ffan(mnew,gama) ;
                     deriv = (rnew-rold)/(mnew-mold) ;
                     rold = rnew ;
                     mold = mnew ;
                     mnew = mold + (rstar - rold)/deriv ;
                     if (mnew < 0.0) mnew = .05 ;
               }
               if (inopt1 == 2) mach1 = mold ;
               else mach2 = mold ;
               break;
       }
       case 5: {                /* pt2/pt1 given */
               if (ratio > 1.0) {
                   erflo =1 ;
                   fl_set_input(outf,"pt2>pt1");
                   return ;
               }
               if(inopt1 == 2) ratio = 1.0/ratio ;
               rstar = get_ptfan(mread,gama)*ratio ;
               if (rstar <= 1.0) {
                   erflo =1 ;
                   fl_set_input(outf,"Choked");
                   return ;
               }
               rold = 1.0 ;
               mold = 1.0 ;
               mnew = mread ;
               iter = 0 ;
               while (fabs(rstar - rold) > .0001) {
                     ++ iter ;
                     rnew = get_ptfan(mnew,gama) ;
                     deriv = (rnew-rold)/(mnew-mold) ; 
                     rold = rnew ;
                     mold = mnew ;
                     mnew = mold + (rstar - rold)/deriv ;
                     if (mnew < 0.0) mnew = .05 ;
               }
               if (inopt1 == 2) mach1 = mold ;
               else mach2 = mold ;
               break;
       }
       case 6: {                /* V2/V1 given */
               if (mread < 1.0 && ratio < 1.0) {
                   erflo =1 ;
                   fl_set_input(outf,"V2<V1");
                   return ;
               }
               if (mread > 1.0 && ratio > 1.0) {
                   erflo =1 ;
                   fl_set_input(outf,"V2>V1");
                   return ;
               }
               if(inopt1 == 2) ratio = 1.0/ratio ;
               rstar = get_vfan(mread,gama)*ratio ;
               if (mread < 1.0 && rstar > 1.0) {
                   erflo =1 ;
                   fl_set_input(outf,"Choked");
                   return ;
               }
               if (mread > 1.0 && rstar < 1.0) {
                   erflo =1 ;
                   fl_set_input(outf,"Choked");
                   return ;
               }
               machs = 2.0*rstar*rstar/(gp1-gm1*rstar*rstar) ;
               if (inopt1 == 2)  mach1 = sqrt(machs) ;
               else mach2 = sqrt(machs) ;
               break;
       }
       case 7: {                /* rho2/rho1 given */
               if (mread < 1.0 && ratio > 1.0) {
                   erflo =1 ;
                   fl_set_input(outf,"r2>r1");
                   return ;
               }
               if (mread > 1.0 && ratio < 1.0) {
                   erflo =1 ;
                   fl_set_input(outf,"r2<r1");
                   return ;
               }
               if(inopt1 == 2) ratio = 1.0/ratio ;
               ratio = 1.0/ratio ;
               rstar = get_vfan(mread,gama)*ratio ;
               if (mread < 1.0 && rstar > 1.0) {
                   erflo =1 ;
                   fl_set_input(outf,"Choked");
                   return ;
               }
               if (mread > 1.0 && rstar < 1.0) {
                   erflo =1 ;
                   fl_set_input(outf,"Choked");
                   return ;
               }
               machs = 2.0*rstar*rstar/(gp1-gm1*rstar*rstar) ;
               if (inopt1 == 2)  mach1 = sqrt(machs) ;
               else mach2 = sqrt(machs) ;
               break;
       }
       case 8: {                /* length given */
               if(inopt1 == 1) rstar = get_lfan(mread,gama) - ratio ;
               if(inopt1 == 2) rstar = ratio + get_lfan(mread,gama) ;
               if (rstar < 0.0 ) {
                   erflo =1 ;
                   fl_set_input(outf,"Choked");
                   return ;
               }
               if (mread > 1.0 && rstar > .8) {
                   erflo =1 ;
                   fl_set_input(outf,"Impossible");
                   return ;
               }
               rold = 0.0 ;
               mold = 1.0 ;
               mnew = mread ;
               while (fabs(rstar - rold) > .0001) {
                     rnew = get_lfan(mnew,gama) ;
                     deriv = (rnew-rold)/(mnew-mold) ;
                     rold = rnew ;
                     mold = mnew ;
                     mnew = mold + (rstar - rold)/deriv ;
                     if (mnew < 0.0) mnew = .05 ;
               }
               if (inopt1 == 2) mach1 = mold ;
               else mach2 = mold ;
               break;
       }
   }

   fanno(mach1,mach2,gama) ;  
   return ;
}

float get_pfan(float machin, float gama)
/* Utility to get the fanno pressure ratio given the mach number */
{
    float number,gp1,gm1;
    gp1 = gama + 1.0 ;
    gm1 = gama - 1.0 ;
    number = sqrt(gp1/(2.0+gm1*machin*machin))/machin ;

    return(number) ;
}

float get_vfan(float machin, float gama)
/* Utility to get the fanno velocity ratio given the mach number */
{
    float number,gp1,gm1;
    gp1 = gama + 1.0 ;
    gm1 = gama - 1.0 ;
    number = machin * sqrt(gp1/(2.0+gm1*machin*machin)) ;

    return(number) ;
}

float get_tfan(float machin, float gama)
/* Utility to get the fanno temperature ratio given the mach number */
{
    float number,gp1,gm1;
    gp1 = gama + 1.0 ;
    gm1 = gama - 1.0 ;
    number = gp1/(2.0+gm1*machin*machin) ;

    return(number) ;
}

float get_ptfan(float machin, float gama)
/* Utility to get the fanno total pressure ratio given the mach number */
{
    float number,fac1,gp1,gm1;
    gp1 = gama + 1.0 ;
    gm1 = gama - 1.0 ;
    fac1 = (2.0+gm1*machin*machin)/gp1 ;
    number = sqrt(pow(fac1,gp1/gm1))/machin ;

    return(number) ;
}

float get_ffan(float machin, float gama)
/* Utility to get the fanno impulse function ratio given the mach number */
{
    float number,fac1,gp1,gm1;
    gp1 = gama + 1.0 ;
    gm1 = gama - 1.0 ;
    fac1 = gp1*(2.0+gm1*machin*machin) ;
    number = (1.0+gama*machin*machin)/(machin*sqrt(fac1)) ;

    return(number) ;
}

float get_lfan(float machin, float gama)
/* Utility to get the fanno l/D given the mach number */
{
    float number,fac1,fac2,gp1,gm1;
    gp1 = gama + 1.0 ;
    gm1 = gama - 1.0 ;
    fac1 = gp1*machin*machin/(2.0+gm1*machin*machin) ;
    fac2 = (1.0-machin*machin)/(gama*machin*machin);
    number = fac2 + .5*gp1/gama*log(fac1) ;

    return(number) ;
}

void fanno (float machin, float machout, float gama)
{                           /* get fanno flow ratios */
   float get_pfan(float,float);
   float get_vfan(float,float);
   float get_tfan(float,float);
   float get_ptfan(float,float);
   float get_ffan(float,float);
   float get_lfan(float,float);
 
   prat  = get_pfan(machout,gama) /get_pfan(machin,gama) ;
   vrat  = get_vfan(machout,gama) /get_vfan(machin,gama) ;
   rhorat = 1.0/vrat ;
   trat  = get_tfan(machout,gama) /get_tfan(machin,gama) ;
   ptrat = get_ptfan(machout,gama)/get_ptfan(machin,gama) ;
   frat  = get_ffan(machout,gama)/get_ffan(machin,gama) ;
   leng  = get_lfan(machin,gama) - get_lfan(machout,gama) ;
}
   
/*      ------------    Atmospheric Calculator  ------- */

void atm_swtch (FL_OBJECT *obj, long arg)
/*   switching routine for Atmosphere tables - uses radio button*/
{
  void anlya () ;

  gama = 1.4 ;          /* get gamma */
  anlya () ;
  switch (arg)
  {
                      /* output */
   case 1: fl_putr (outa,ps,2) ; break ;    /* static press */
   case 2: fl_putr (outa,ptrat,0) ; break ; /* total press  */
   case 3: fl_putr (outa,ts,2) ; break ;    /* static temp  */
   case 4: fl_putr (outa,ttrat,2) ; break ; /* total temp   */
   case 5: fl_putr (outa,q0,0) ; break ;    /* q  */
   case 6: fl_putr (outa,u0,0) ; break ;    /* uvel  */
   case 7: fl_putr (outa,mu,4) ; break ;    /* viscosity  */
   case 8: fl_putr (outa,rey,0) ; break ;   /* reynolds number */
  }
}

void create_aform()
{
  FL_OBJECT *obj;

  aform = fl_bgn_form(FL_UP_BOX,190.0,250.0);

               /* Input */
  obj = fl_add_box(FL_NO_BOX,5.0,220.0,90.0,25.0,"Mach");
  ina1 = fl_add_input(FL_NORMAL_INPUT,5.,190.,90.,25.,"") ;
     fl_set_object_color(ina1,7,7);
     fl_set_object_lsize(ina1,FL_MEDIUM_SIZE);
  obj = fl_add_box(FL_NO_BOX,100.0,220.0,90.0,25.0,"Alt -ft");
  ina2 = fl_add_input(FL_NORMAL_INPUT,100.,190.,85.,25.,"") ;
     fl_set_object_color(ina2,7,7);
     fl_set_object_lsize(ina2,FL_MEDIUM_SIZE);

             /* Choice of Outputs */
  obj = fl_add_box(FL_FLAT_BOX,5.0,30.0,180.0,155.0,"");
    fl_set_object_color(obj,7,7);
    obj = fl_add_lightbutton(FL_RADIO_BUTTON,5.0,155.0,90.0,25.0,"p-lb/sft");
       fl_set_object_callback(obj,atm_swtch,1);
    obj = fl_add_lightbutton(FL_RADIO_BUTTON,100.0,155.0,90.0,25.0,"pt-lb/sft");
       fl_set_object_callback(obj,atm_swtch,2);
    obj = fl_add_lightbutton(FL_RADIO_BUTTON,5.0,125.0,90.0,25.0,"Ts - R");
       fl_set_object_callback(obj,atm_swtch,3);
    obj = fl_add_lightbutton(FL_RADIO_BUTTON,100.0,125.0,90.0,25.0,"Tt - R");
       fl_set_object_callback(obj,atm_swtch,4);
    obj = fl_add_lightbutton(FL_RADIO_BUTTON,5.0,95.0,90.0,25.0,"q-lb/sft");
       fl_set_object_callback(obj,atm_swtch,5);
    obj = fl_add_lightbutton(FL_RADIO_BUTTON,100.0,95.0,90.0,25.0,"u-ft/s");
       fl_set_object_callback(obj,atm_swtch,6);
    obj = fl_add_lightbutton(FL_RADIO_BUTTON,5.0,65.0,90.0,25.0,"nu-sft/s");
       fl_set_object_callback(obj,atm_swtch,7);
    obj = fl_add_lightbutton(FL_RADIO_BUTTON,100.0,65.0,90.0,25.0,"Re/ft");
       fl_set_object_callback(obj,atm_swtch,8);
             /* Output window */
  outa = fl_add_input(FL_NORMAL_INPUT,45.,35.,100.,25.,"") ;
       fl_set_input_color(outa,0,7);
       fl_set_object_lsize(outa,FL_MEDIUM_SIZE);

              /* Exit Button */
  obj = fl_add_button(FL_NORMAL_BUTTON,45.0,5.0,90.0,25.0,"Exit");
         fl_set_object_callback(obj,clos_swtch,8); 
         fl_set_object_color(obj,1,1);
         fl_set_object_lcol(obj,7);

  fl_end_form();
}

void anlya ()                           /* analysis for atmosphere */
{
    float get_tisen(float,float) ;
    float get_pisen(float,float) ;

    float rgas = 1718.;
    if (gama <1.01) gama = 1.01 ;
    gm1 = gama - 1.0 ;
    gp1 = gama + 1.0;
    mach1 = fl_getr(ina1) ;
    alt   = fl_getr(ina2);
    if (alt <= 36000. ) {
           ts = 518.6 - 3.56 * alt / 1000. ;
           ps = 2116. * pow(ts/518.6, 5.256) ;
    }
    if (alt <= 82300. && alt > 36000.) {
           ts = 389.98 ;
           ps = 2116. * .2236 * exp((36000-alt)/(53.35*389.98)) ;
    }
    if (alt < 160000. && alt > 82300.) {
           ts = 389.9 + 1.645 * (alt-82300.) / 1000. ;
           ps = 2116. * .0245 * exp((82300-alt)/(53.35*389.98)) ;
    }
    u0 = mach1 * sqrt(gama*rgas*ts) ;
    q0 = gama/2.0*mach1*mach1*ps ;
    ttrat = ts / get_tisen(mach1,gama) ;
    ptrat = ps / get_pisen(mach1,gama) ;
                        /* viscosity from Sutherland's law */
    mu = .0020886*.1716*690.6/(ts + 199.)*pow(ts/491.6,1.5) ;
    mu = mu * rgas *ts /ps ;
    rey =  u0 * 1000. / mu ;
    return ;
}
 
void init_forms ()
{                                       /* initial value for variables */
   fl_putr (gami,1.4,2) ;
   fl_putr (ini1,2.0,2) ;
   fl_putr (gamn,1.4,2) ;
   fl_putr (inn1,2.0,2) ;
   fl_putr (gamo,1.4,2) ;
   fl_putr (ino1,2.0,2) ;
   fl_putr (ino2,10.0,2) ;
   fl_putr (game,1.4,2) ;
   fl_putr (ine1,2.0,2) ;
   fl_putr (ine2,10.0,2) ;
   fl_putr (gamr,1.4,2) ;
   fl_putr (inr1,.5,2) ;
   fl_putr (inr2,.6,2) ;
   fl_putr (gamf,1.4,2) ;
   fl_putr (inf1,.5,2) ;
   fl_putr (inf2,.6,2) ;
   fl_putr (ina1,2.0,2) ;
   fl_putr (ina2,15000.,0) ;
}

main(int argc, char *argv[])
{
  FL_OBJECT *obj;
  long dev;
  short val;

                     /* initialize the forms */
   fl_initialize(argv[0], "FormDemo", 0, 0 ,&argc, argv);
   create_mainform ();
   create_hlpform ();
   create_diagform ();
   create_nform ();
   create_iform ();
   create_oform ();
   create_eform ();
   create_rform ();
   create_fform ();
   create_aform ();
   init_forms () ;
   fl_show_form(mainform,FL_PLACE_ASPECT,FL_FULLBORDER,"AERO");
  fl_set_atclose(at_close,"") ;

  while (1)                      /* forms loop */
  {
    obj = fl_do_forms();       
    if (obj == killbut) exit(0);
  }
}
