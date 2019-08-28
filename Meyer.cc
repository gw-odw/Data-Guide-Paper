// Wavelet Analysis Tool
//--------------------------------------------------------------------
// Implementation of
// Meyer wavelets using Fast Wavelet Transform
// References:
//--------------------------------------------------------------------
//$Id: Meyer.hh,v 0.2 2001/08/06 19:37:00 klimenko Exp $

#define MEYER_CC

#include "Meyer.hh"
#include <iostream>
#include <cstdio>

//namespace datacondAPI {
//namespace wat {

extern const double mey[16]=
{ 7.4458559231880628e-01,   4.4459300275757724e-01,  -3.5087555656258346e-02,
 -1.3284520043622938e-01,   3.0655091960824263e-02,   6.3739024322801596e-02,
 -2.4348745906078023e-02,  -3.2130793990211758e-02,   1.7423434103729693e-02,
  1.5270015130934803e-02,  -1.1061496392513451e-02,  -6.3877183184971563e-03,
  6.0458140973233040e-03,   2.2025341009110021e-03,  -2.7046721246437250e-03,
 -6.0115023435160925e-04};

extern const double mey004[32]=
{  7.43750353e-01,   4.44094660e-01,  -3.50482247e-02,  -1.32696489e-01,
   3.06212119e-02,   6.36672446e-02,  -2.43217736e-02,  -3.20940567e-02,
   1.74038765e-02,   1.52509309e-02,  -1.10446159e-02,  -6.38673997e-03,
   6.04547967e-03,   2.19476459e-03,  -2.70215388e-03,  -5.78163043e-04,
   8.59511073e-04,   1.61178251e-04,  -9.35363005e-05,  -1.39924765e-04,
  -7.55811503e-05,   1.49736854e-04,   2.44738777e-05,  -1.02790685e-04,
   3.73859490e-05,   3.23295640e-05,  -4.64169099e-05,   1.50170027e-05,
   2.01182002e-05,  -2.44497410e-05,   6.16393792e-06,   1.04921553e-05};

extern const double mey008[48]=
{  7.33609639e-01,   4.46956642e-01,  -2.58537823e-02,  -1.40705138e-01,
   2.39953217e-02,   7.52666880e-02,  -2.11749933e-02,  -4.52098615e-02,
   1.77478786e-02,   2.78607152e-02,  -1.41068440e-02,  -1.69993874e-02,
   1.06121463e-02,   1.00922602e-02,  -7.53699278e-03,  -5.78259209e-03,
   5.03950041e-03,   3.19710122e-03,  -3.16307359e-03,  -1.72272472e-03,
   1.85945070e-03,   9.25826653e-04,  -1.02391321e-03,  -5.12999856e-04,
   5.31379384e-04,   3.00631034e-04,  -2.64764748e-04,  -1.85201368e-04,
   1.31470288e-04,   1.15431643e-04,  -6.82864331e-05,  -6.94283079e-05,
   3.79892424e-05,   3.88121845e-05,  -2.19363081e-05,  -1.97440301e-05,
   1.23471763e-05,   9.18041839e-06,  -6.45040920e-06,  -4.13481950e-06,
   3.11352396e-06,   2.05800311e-06,  -1.49196971e-06,  -1.22615994e-06,
   8.12801410e-07,   7.67372574e-07,  -5.15136932e-07,  -4.06233507e-07};

extern const double mey016[64]=
{  7.26060766e-01,   4.48510741e-01,  -1.87063583e-02,  -1.45182691e-01,
   1.79824715e-02,   8.21488422e-02,  -1.68369560e-02,  -5.37407521e-02,
   1.53534264e-02,   3.71819869e-02,  -1.36345964e-02,  -2.62900015e-02,
   1.17907408e-02,   1.86830184e-02,  -9.92827064e-03,  -1.32245702e-02,
   8.14003237e-03,   9.27736602e-03,  -6.49849436e-03,  -6.43454314e-03,
   5.05236309e-03,   4.41002788e-03,  -3.82653735e-03,  -2.99022990e-03,
   2.82480333e-03,   2.01133358e-03,  -2.03436735e-03,  -1.34748651e-03,
   1.43122726e-03,   9.03470490e-04,  -9.85466744e-04,  -6.09099247e-04,
   6.65769880e-04,   4.14298193e-04,  -4.42743822e-04,  -2.84604394e-04,
   2.90941691e-04,   1.97138048e-04,  -1.89729736e-04,  -1.37177063e-04,
   1.23292062e-04,   9.54358756e-05,  -8.01004803e-05,  -6.60692637e-05,
   5.21278478e-05,   4.53290749e-05,  -3.40029351e-05,  -3.07304146e-05,
   2.22340207e-05,   2.05603112e-05,  -1.45771375e-05,  -1.35893549e-05,
   9.58175075e-06,   8.90067120e-06,  -6.30288149e-06,  -5.79549494e-06,
   4.13268817e-06,   3.75088503e-06,  -2.69342158e-06,  -2.40111942e-06,
   1.75266393e-06,   1.51368900e-06,  -1.15405288e-06,  -9.47428710e-07};

extern const double mey032[96]=
{  7.20585503e-01,   4.49322176e-01,  -1.33876233e-02,  -1.47564003e-01,
   1.31180419e-02,   8.59468746e-02,  -1.26808775e-02,  -5.87187863e-02,
   1.20935038e-02,   4.30446119e-02,  -1.13786310e-02,  -3.27126251e-02,
   1.05628353e-02,   2.53419265e-02,  -9.67489952e-03,  -1.98236902e-02,
   8.74412541e-03,   1.55687636e-02,  -7.79876613e-03,  -1.22309446e-02,
   6.86468095e-03,   9.58904815e-03,  -5.96427568e-03,  -7.49122252e-03,
   5.11576757e-03,   5.82666639e-03,  -4.33280026e-03,  -4.51048065e-03,
   3.62441235e-03,   3.47522836e-03,  -2.99532338e-03,  -2.66608068e-03,
   2.44645878e-03,   2.03792774e-03,  -1.97561400e-03,  -1.55355878e-03,
   1.57816957e-03,   1.18240002e-03,  -1.24780015e-03,  -8.99529934e-04,
   9.77145013e-04,   6.84842091e-04,  -7.58411779e-04,  -5.22310303e-04,
   5.83878332e-04,   3.99340552e-04,  -4.46261029e-04,  -3.06191991e-04,
   3.38941838e-04,   2.35446834e-04,  -2.56079543e-04,  -1.81522688e-04,
   1.92647355e-04,   1.40242548e-04,  -1.44429583e-04,  -1.08485956e-04,
   1.07985946e-04,   8.39294256e-05,  -8.05779149e-05,  -6.48578157e-05,
   6.00598599e-05,   5.00147104e-05,  -4.47575881e-05,  -3.84709859e-05,
   3.33650554e-05,   2.95145229e-05,  -2.48751954e-05,  -2.25766867e-05,
   1.85348434e-05,   1.72006556e-05,  -1.38002625e-05,  -1.30344747e-05,
   1.02805058e-05,   9.82123371e-06,  -7.67887314e-06,  -7.37168717e-06,
   5.75432885e-06,   5.52911299e-06,  -4.31302644e-06,  -4.14866747e-06,
   3.21670349e-06,   3.10218711e-06,  -2.38316834e-06,  -2.29613316e-06,
   1.76612899e-06,   1.67878805e-06,  -1.32540443e-06,  -1.22471389e-06,
   1.01052828e-06,   9.08116537e-07,  -7.69114071e-07,  -6.88574605e-07,
   5.67889970e-07,   5.20825336e-07,  -4.02604215e-07,  -3.76740644e-07};

extern const double mey064[128]=
{  7.16664739e-01,   4.49737013e-01,  -9.52509311e-03,  -1.48794259e-01,
   9.42718353e-03,   8.79506343e-02,  -9.26626565e-03,  -6.14289433e-02,
   9.04568038e-03,   4.63727568e-02,  -8.76996656e-03,  -3.65540285e-02,
   8.44469228e-03,   2.95812182e-02,  -8.07624332e-03,  -2.43407973e-02,
   7.67160099e-03,   2.02446765e-02,  -7.23813107e-03,  -1.69529251e-02,
   6.78338198e-03,   1.42549718e-02,  -6.31487811e-03,  -1.20128486e-02,
   5.83990528e-03,   1.01317354e-02,  -5.36530749e-03,  -8.54364745e-03,
   4.89732421e-03,   7.19792811e-03,  -4.44148387e-03,  -6.05549731e-03,
   4.00254383e-03,   5.08527746e-03,  -3.58445343e-03,  -4.26192058e-03,
   3.19032657e-03,   3.56431865e-03,  -2.82243293e-03,  -2.97459250e-03,
   2.48222706e-03,   2.47739105e-03,  -2.17042247e-03,  -2.05941104e-03,
   1.88709354e-03,   1.70907625e-03,  -1.63177647e-03,  -1.41631659e-03,
   1.40355229e-03,   1.17239254e-03,  -1.20111920e-03,  -9.69732697e-04,
   1.02287367e-03,   8.01783161e-04,  -8.67009056e-04,  -6.62884201e-04,
   7.31617418e-04,   5.48181040e-04,  -6.14769429e-04,  -4.53554313e-04,
   5.14559732e-04,   3.75545294e-04,  -4.29129354e-04,  -3.11263665e-04,
   3.56689096e-04,   2.58289569e-04,  -2.95556490e-04,  -2.14593685e-04,
   2.44195619e-04,   1.78487296e-04,  -2.01237488e-04,  -1.48590646e-04,
   1.65470505e-04,   1.23795972e-04,  -1.35814182e-04,  -1.03213263e-04,
   1.11300869e-04,   8.61101801e-05,  -9.10784432e-05,  -7.18692981e-05,
   7.44231322e-05,   5.99739223e-05,  -6.07397066e-05,  -5.00101480e-05,
   4.95379001e-05,   4.16610140e-05,  -4.03973174e-05,  -3.46804093e-05,
   3.29447468e-05,   2.88579911e-05,  -2.68558874e-05,  -2.39982035e-05,
   2.18698321e-05,   1.99247507e-05,  -1.77927861e-05,  -1.64983720e-05,
   1.44792262e-05,   1.36240347e-05,  -1.18022353e-05,  -1.12355047e-05,
   9.63651184e-06,   9.26887294e-06,  -7.86574972e-06,  -7.64844191e-06,
   6.40251358e-06,   6.29662555e-06,  -5.19694529e-06,  -5.15597698e-06,
   4.22242785e-06,   4.19969775e-06,  -3.44991009e-06,  -3.41878746e-06,
   2.83439116e-06,   2.79757727e-06,  -2.32528820e-06,  -2.30118524e-06,
   1.88884962e-06,   1.88665321e-06,  -1.51899775e-06,  -1.52596025e-06,
   1.22477622e-06,   1.21732860e-06,  -1.00615445e-06,  -9.73026670e-07,
   8.41733183e-07,   7.95446067e-07,  -7.00115852e-07,  -6.65015290e-07};

extern const double mey128[192]=
{  7.13874857e-01,   4.49946782e-01,  -6.75634361e-03,  -1.49419865e-01,
   6.72126119e-03,   8.89810886e-02,  -6.66318178e-03,  -6.28463092e-02,
   6.58271363e-03,   4.81526973e-02,  -6.48072078e-03,  -3.86665398e-02,
   6.35829589e-03,   3.19915158e-02,  -6.21671757e-03,  -2.70103180e-02,
   6.05741630e-03,   2.31321538e-02,  -5.88196124e-03,  -2.00154955e-02,
   5.69205682e-03,   1.74492695e-02,  -5.48952611e-03,  -1.52960492e-02,
   5.27626983e-03,   1.34625138e-02,  -5.05421361e-03,  -1.18830206e-02,
   4.82526746e-03,   1.05099617e-02,  -4.59131003e-03,  -9.30786174e-03,
   4.35418637e-03,   8.24964702e-03,  -4.11569602e-03,  -7.31421627e-03,
   3.87756010e-03,   6.48480242e-03,  -3.64137952e-03,  -5.74782559e-03,
   3.40860795e-03,   5.09207359e-03,  -3.18055166e-03,  -4.50812305e-03,
   2.95838430e-03,   3.98794154e-03,  -2.74315314e-03,  -3.52461175e-03,
   2.53576474e-03,   3.11212215e-03,  -2.33696174e-03,  -2.74519036e-03,
   2.14731404e-03,   2.41911564e-03,  -1.96723589e-03,  -2.12967297e-03,
   1.79701675e-03,   1.87305269e-03,  -1.63684203e-03,  -1.64582811e-03,
   1.48679163e-03,   1.44492316e-03,  -1.34682759e-03,  -1.26756542e-03,
   1.21679416e-03,   1.11123371e-03,  -1.09644184e-03,  -9.73622283e-04,
   9.85463097e-04,   8.52632074e-04,  -8.83516198e-04,  -7.46376303e-04,
   7.90224945e-04,   6.53176176e-04,  -7.05166027e-04,  -5.71534445e-04,
   6.27867363e-04,   5.00098286e-04,  -5.57829147e-04,  -4.37634821e-04,
   4.94555708e-04,   3.83030948e-04,  -4.37574603e-04,  -3.35305600e-04,
   3.86431130e-04,   2.93610816e-04,  -3.40670064e-04,  -2.57209815e-04,
   2.99828230e-04,   2.25443881e-04,  -2.63449743e-04,  -1.97711592e-04,
   2.31112201e-04,   1.73472226e-04,  -2.02440347e-04,  -1.52261530e-04,
   1.77095480e-04,   1.33696306e-04,  -1.54752496e-04,  -1.17456008e-04,
   1.35088212e-04,   1.03253143e-04,  -1.17792828e-04,  -9.08160149e-05,
   1.02592849e-04,   7.98955829e-05,  -8.92619378e-05,  -7.02846186e-05,
   7.76079954e-05,   6.18255623e-05,  -6.74483122e-05,  -5.43952605e-05,
   5.85964038e-05,   4.78783358e-05,  -5.08723664e-05,  -4.21527247e-05,
   4.41250033e-05,   3.70991295e-05,  -3.82421809e-05,  -3.26225654e-05,
   3.31376531e-05,   2.86623967e-05,  -2.87261562e-05,  -2.51790476e-05,
   2.49103579e-05,   2.21291417e-05,  -2.15914552e-05,  -1.94526151e-05,
   1.86916430e-05,   1.70835602e-05,  -1.61648846e-05,  -1.49729930e-05,
   1.39841985e-05,   1.30999537e-05,  -1.21172445e-05,  -1.14591373e-05,
   1.05137757e-05,   1.00368267e-05,  -9.11673804e-06,  -8.79868318e-06,
   7.88522928e-06,   7.70116876e-06,  -6.80573957e-06,  -6.71480408e-06,
   5.87988418e-06,   5.83568375e-06,  -5.10040617e-06,  -5.07345765e-06,
   4.43901459e-06,   4.42755830e-06,  -3.85783858e-06,  -3.87523967e-06,
   3.33270839e-06,   3.38320811e-06,  -2.86468813e-06,  -2.93105633e-06,
   2.46807056e-06,   2.52292752e-06,  -2.14661427e-06,  -2.17562281e-06,
   1.88159042e-06,   1.89493599e-06,  -1.64342181e-06,  -1.66378420e-06,
   1.41512696e-06,   1.45391789e-06,  -1.20399650e-06,  -1.24942327e-06,
   1.02971383e-06,   1.05844635e-06,  -9.00702688e-07,  -9.01352903e-07,
   8.02270160e-07,   7.87109450e-07,  -7.08329725e-07,  -7.01454992e-07,
   6.04918371e-07,   6.18647449e-07,  -5.01936480e-07,  -5.24998941e-07,
   4.21323580e-07,   4.30629196e-07,  -3.73453754e-07,  -3.57652076e-07,
   3.45320152e-07,   3.16580472e-07,  -3.12293890e-07,  -2.94519438e-07,
   2.61672192e-07,   2.66932032e-07,  -2.04445070e-07,  -2.21192051e-07};

extern const double mey256[256]=
{  7.11895915e-01,   4.50052264e-01,  -4.78497187e-03,  -1.49735366e-01,
   4.77249670e-03,   8.95037939e-02,  -4.75175254e-03,  -6.35715771e-02,
   4.72283975e-03,   4.90741318e-02,  -4.68592571e-03,  -3.97760905e-02,
   4.64123155e-03,   3.32796024e-02,  -4.58900682e-03,  -2.84659765e-02,
   4.52951548e-03,   2.47432011e-02,  -4.46304525e-03,  -2.17687091e-02,
   4.38992835e-03,   1.93305783e-02,  -4.31055029e-03,  -1.72907313e-02,
   4.22533483e-03,   1.55553974e-02,  -4.13471695e-03,  -1.40586841e-02,
   4.03912748e-03,   1.27529294e-02,  -3.93900127e-03,  -1.16027913e-02,
   3.83479696e-03,   1.05815055e-02,  -3.72700512e-03,  -9.66844238e-03,
   3.61613270e-03,   8.84745277e-03,  -3.50267588e-03,  -8.10570245e-03,
   3.38710476e-03,   7.43283362e-03,  -3.26987177e-03,  -6.82036590e-03,
   3.15143203e-03,   6.26127812e-03,  -3.03225220e-03,  -5.74971233e-03,
   2.91279589e-03,   5.28074500e-03,  -2.79349761e-03,  -4.85019185e-03,
   2.67474873e-03,   4.45444318e-03,  -2.55690733e-03,  -4.09034238e-03,
   2.44032007e-03,   3.75511180e-03,  -2.32533262e-03,  -3.44630841e-03,
   2.21227676e-03,   3.16178175e-03,  -2.10144605e-03,  -2.89961896e-03,
   1.99308348e-03,   2.65808659e-03,  -1.88739305e-03,  -2.43559083e-03,
   1.78456332e-03,   2.23066657e-03,  -1.68477945e-03,  -2.04198243e-03,
   1.58821189e-03,   1.86833741e-03,  -1.49499352e-03,  -1.70863655e-03,
   1.40520867e-03,   1.56185713e-03,  -1.31890598e-03,  -1.42702822e-03,
   1.23612311e-03,   1.30323532e-03,  -1.15689982e-03,  -1.18963786e-03,
   1.08126754e-03,   1.08547586e-03,  -1.00922727e-03,  -9.90053688e-04,
   9.40739287e-04,   9.02712742e-04,  -8.75736458e-04,  -8.22816264e-04,
   8.14149362e-04,   7.49758197e-04,  -7.55919615e-04,  -6.82984093e-04,
   7.00989599e-04,   6.22000499e-04,  -6.49280388e-04,  -5.66360968e-04,
   6.00681390e-04,   5.15640447e-04,  -5.55063515e-04,  -4.69421588e-04,
   5.12304059e-04,   4.27304738e-04,  -4.72299723e-04,  -3.88929801e-04,
   4.34955986e-04,   3.53986405e-04,  -4.00164607e-04,  -3.22200558e-04,
   3.67792818e-04,   2.93309586e-04,  -3.37696007e-04,  -2.67048920e-04,
   3.09742089e-04,   2.43162506e-04,  -2.83824007e-04,  -2.21425047e-04,
   2.59848581e-04,   2.01652528e-04,  -2.37713474e-04,  -1.83689204e-04,
   2.17295879e-04,   1.67382873e-04,  -1.98464686e-04,  -1.52571976e-04,
   1.81104372e-04,   1.39096329e-04,  -1.65127035e-04,  -1.26819686e-04,
   1.50460798e-04,   1.15640574e-04,  -1.37026367e-04,  -1.05479612e-04,
   1.24725331e-04,   9.62550897e-05,  -1.13451977e-04,  -8.78703983e-05,
   1.03116857e-04,   8.02250687e-05,  -9.36585297e-05,  -7.32376525e-05,
   8.50317076e-05,   6.68568635e-05,  -7.71835855e-05,  -6.10491963e-05,
   7.00419362e-05,   5.57748051e-05,  -6.35267576e-05,  -5.09752111e-05,
   5.75736906e-05,   4.65846222e-05,  -5.21456430e-05,  -4.25530751e-05,
   4.72208342e-05,   3.88578298e-05,  -4.27690495e-05,  -3.54912276e-05,
   3.87396764e-05,   3.24367969e-05,  -3.50733108e-05,  -2.96571742e-05,
   3.17251496e-05,   2.71056237e-05,  -2.86766020e-05,  -2.47493665e-05,
   2.59233336e-05,   2.25811494e-05,  -2.34515310e-05,  -2.06072631e-05,
   2.12259555e-05,   1.88237971e-05,  -1.92015749e-05,  -1.72046971e-05,
   1.73469870e-05,   1.57134092e-05,  -1.56560663e-05,  -1.43263233e-05,
   1.41360474e-05,   1.30444452e-05,  -1.27838306e-05,  -1.18815095e-05,
   1.15740798e-05,   1.08403188e-05,  -1.04708976e-05,  -9.90087789e-06,
   9.45129366e-06,   9.03210675e-06,  -8.51687574e-06,  -8.21534571e-06,
   7.68197656e-06,   7.45608240e-06,  -6.95000216e-06,  -6.77211429e-06,
   6.30157083e-06,   6.16993176e-06,  -5.70622823e-06,  -5.63289129e-06,
   5.14595364e-06,   5.13296307e-06,  -4.62688686e-06,  -4.65426709e-06,
   4.16748961e-06,   4.20482738e-06,  -3.77492600e-06,  -3.80475880e-06,
   3.43323180e-06,   3.46267165e-06,  -3.11505790e-06,  -3.16386391e-06,
   2.80520299e-06,   2.88208567e-06,  -2.51236489e-06,  -2.60309036e-06,
   2.25732477e-06,   2.33640219e-06,  -2.04934863e-06,  -2.10351500e-06,
   1.87437601e-06,   1.91430744e-06,  -1.70678134e-06,  -1.75524474e-06,
   1.53292300e-06,   1.60115165e-06,  -1.36290960e-06,  -1.43877127e-06,
   1.21879802e-06,   1.27853951e-06,  -1.11100740e-06,  -1.14279008e-06,
   1.02651944e-06,   1.04217548e-06,  -9.40650171e-07,  -9.63873785e-07,
   8.40608522e-07,   8.83366419e-07,  -7.37271313e-07,  -7.88001215e-07,
   6.53389032e-07,   6.88770812e-07,  -6.00007052e-07,  -6.08521347e-07,
   5.64672395e-07,   5.58376709e-07,  -5.23211431e-07,  -5.25948504e-07,
   4.63293760e-07,   4.87116590e-07,  -3.96211976e-07,  -4.29594843e-07,
   3.45091793e-07,   3.64711991e-07,  -3.21317204e-07,  -3.15622573e-07};

extern const double mey384[384]=
{  7.11017584e-01,   4.50087409e-01,  -3.90864520e-03,  -1.49840731e-01,
   3.90206191e-03,   8.96790649e-02,  -3.89084227e-03,  -6.38159860e-02,
   3.87490944e-03,   4.93863601e-02,  -3.85443835e-03,  -4.01543853e-02,
   3.82973748e-03,   3.37220018e-02,  -3.80101268e-03,  -2.89704313e-02,
   3.76824883e-03,   2.53074870e-02,  -3.73132705e-03,  -2.23902353e-02,
   3.69025970e-03,   2.00063033e-02,  -3.64530741e-03,  -1.80172933e-02,
   3.59686042e-03,   1.63293518e-02,  -3.54520195e-03,  -1.48766360e-02,
   3.49038937e-03,   1.36114551e-02,  -3.43237106e-03,  -1.24982537e-02,
   3.37122106e-03,   1.15099762e-02,  -3.30725584e-03,  -1.06258357e-02,
   3.24091535e-03,   9.82976528e-03,  -3.17252619e-03,  -9.10914602e-03,
   3.10218270e-03,   8.45375526e-03,  -3.02986363e-03,  -7.85506037e-03,
   2.95566681e-03,   7.30590531e-03,  -2.87992593e-03,  -6.80042517e-03,
   2.80309157e-03,   6.33392168e-03,  -2.72549458e-03,  -5.90256016e-03,
   2.64722721e-03,   5.50299075e-03,  -2.56826006e-03,  -5.13211877e-03,
   2.48867691e-03,   4.78713456e-03,  -2.40879175e-03,  -4.46567920e-03,
   2.32903023e-03,   4.16590632e-03,  -2.24969322e-03,  -3.88631909e-03,
   2.17083836e-03,   3.62549780e-03,  -2.09239740e-03,  -3.38195207e-03,
   2.01441139e-03,   3.15421411e-03,  -1.93714812e-03,  -2.94105418e-03,
   1.86098398e-03,   2.74158159e-03,  -1.78616796e-03,  -2.55511292e-03,
   1.71270364e-03,   2.38092463e-03,  -1.64046686e-03,  -2.21812548e-03,
   1.56944142e-03,   2.06576607e-03,  -1.49983691e-03,  -1.92306771e-03,
   1.43197091e-03,   1.78953434e-03,  -1.36603342e-03,  -1.66482974e-03,
   1.30196917e-03,   1.54853763e-03,  -1.23959572e-03,  -1.44004013e-03,
   1.17883935e-03,   1.33863258e-03,  -1.11985330e-03,  -1.24375656e-03,
   1.06290017e-03,   1.15511541e-03,  -1.00811665e-03,  -1.07255448e-03,
   9.55396053e-04,   9.95823582e-04,  -9.04506526e-04,  -9.24457717e-04,
   8.55327235e-04,   8.57893516e-04,  -8.07966638e-04,  -7.95703775e-04,
   7.62645097e-04,   7.37714244e-04,  -7.19459630e-04,  -6.83884826e-04,
   6.78266519e-04,   6.34073000e-04,  -6.38799615e-04,  -5.87915177e-04,
   6.00906492e-04,   5.44943812e-04,  -5.64666746e-04,  -5.04822423e-04,
   5.30274583e-04,   4.67462815e-04,  -4.97803536e-04,  -4.32906634e-04,
   4.67089022e-04,   4.01089094e-04,  -4.37846574e-04,  -3.71720589e-04,
   4.09907918e-04,   3.44404015e-04,  -3.83339184e-04,  -3.18869979e-04,
   3.58323383e-04,   2.95094182e-04,  -3.34925025e-04,  -2.73179107e-04,
   3.12972551e-04,   2.53117874e-04,  -2.92176448e-04,  -2.34675958e-04,
   2.72365201e-04,   2.17508613e-04,  -2.53603369e-04,  -2.01396171e-04,
   2.36073945e-04,   1.86361512e-04,  -2.19842838e-04,  -1.72551841e-04,
   2.04741179e-04,   1.60002621e-04,  -1.90483310e-04,  -1.48519365e-04,
   1.76902612e-04,   1.37795143e-04,  -1.64069468e-04,  -1.27645963e-04,
   1.52173503e-04,   1.18128308e-04,  -1.41287958e-04,  -1.09420999e-04,
   1.31251892e-04,   1.01589203e-04,  -1.21788072e-04,  -9.44663026e-05,
   1.12738710e-04,   8.77714836e-05,  -1.04183332e-04,  -8.13451855e-05,
   9.63209473e-05,   7.52667260e-05,  -8.92343436e-05,  -6.97362337e-05,
   8.27722214e-05,   6.48387355e-05,  -7.66670221e-05,  -6.04261010e-05,
   7.07706056e-05,   5.62346999e-05,  -6.51720769e-05,  -5.21209285e-05,
   6.00799054e-05,   4.81789030e-05,  -5.55861831e-05,  -4.46224615e-05,
   5.15487245e-05,   4.15493180e-05,  -4.77088656e-05,  -3.88230703e-05,
   4.39271165e-05,   3.61909193e-05,  -4.03009685e-05,  -3.35192549e-05,
   3.70469977e-05,   3.09114062e-05,  -3.42651102e-05,  -2.85896985e-05,
   3.18206326e-05,   2.66596591e-05,  -2.94621038e-05,  -2.49920735e-05,
   2.70569244e-05,   2.33407512e-05,  -2.47091627e-05,  -2.15781546e-05,
   2.26416598e-05,   1.98131916e-05,  -2.09602777e-05,  -1.82733110e-05,
   1.95359940e-05,   1.70687447e-05,  -1.81226984e-05,  -1.60745977e-05,
   1.65928502e-05,   1.50486440e-05,  -1.50552921e-05,  -1.38669855e-05,
   1.37373646e-05,   1.26418718e-05,  -1.27491643e-05,  -1.16038210e-05,
   1.19656476e-05,   1.08658854e-05,  -1.11444373e-05,  -1.03057618e-05,
   1.01614907e-05,   9.68360712e-06,  -9.12892437e-06,  -8.87771523e-06,
   8.27714049e-06,   8.00235344e-06,  -7.71909503e-06,  -7.28989912e-06,
   7.33241173e-06,   6.85511870e-06,  -6.87719924e-06,  -6.57729069e-06,
   6.23172852e-06,   6.21803201e-06,  -5.51026882e-06,  -5.65698534e-06,
   4.94522316e-06,   5.00966549e-06,  -4.65140639e-06,  -4.50960504e-06,
   4.50816730e-06,   4.27264119e-06,  -4.27721386e-06,  -4.17904541e-06,
   3.83829877e-06,   3.99135588e-06,  -3.30706353e-06,  -3.59007062e-06,
   2.91718177e-06,   3.09149959e-06,  -2.78464422e-06,  -2.72991556e-06,
   2.78988692e-06,   2.62184526e-06,  -2.69562226e-06,  -2.64820266e-06,
   2.38253026e-06,   2.57212448e-06,  -1.96710869e-06,  -2.27466668e-06,
   1.68382189e-06,   1.87265988e-06,  -1.64939057e-06,  -1.60086337e-06,
   1.74492415e-06,   1.57625873e-06,  -1.73375546e-06,  -1.68018495e-06,
   1.49713605e-06,   1.67617565e-06,  -1.15209033e-06,  -1.44565738e-06,
   9.33568272e-07,   1.10580727e-06,  -9.58738388e-07,  -8.91709022e-07,
   1.10912250e-06,   9.20648047e-07,  -1.14843279e-06,  -1.07424749e-06,
   9.58269576e-07,   1.11630634e-06,  -6.55978270e-07,  -9.28498858e-07,
   4.76804701e-07,   6.28233782e-07,  -5.38190060e-07,  -4.50811695e-07,
   7.21907128e-07,   5.13721134e-07,  -7.91898561e-07,  -6.98775094e-07,
   6.29976309e-07,   7.69949518e-07,  -3.53680414e-07,  -6.09083564e-07,
   1.98436165e-07,   3.33739875e-07,  -2.81850268e-07,  -1.79363102e-07,
   4.85847782e-07,   2.63576622e-07,  -5.74511016e-07,  -4.68319242e-07,
   4.29779921e-07,   5.57683825e-07,  -1.69312228e-07,  -4.13618087e-07,
   2.86418681e-08,   1.53785781e-07,  -1.25475908e-07,  -1.37261001e-08,
   3.41831748e-07,   1.11150848e-07,  -4.41876227e-07,  -3.28081137e-07,
   3.07626614e-07,   4.28685999e-07,  -5.68117301e-08,  -2.94983550e-07,
  -7.49686950e-08,   4.47029437e-08,  -3.00533793e-08,   8.65554823e-08,
   2.53951480e-07,   1.89772734e-08,  -3.60943678e-07,  -2.43375449e-07,
   2.33093935e-07,   3.50856976e-07,   1.18260092e-08,  -2.23484904e-07,
  -1.38176528e-07,  -2.09699988e-08,   2.81524357e-08,   1.46868453e-07,
   2.00354017e-07,  -3.64049742e-08,  -3.11591993e-07,  -1.92528478e-07,
   1.87653644e-07,   3.04180734e-07,   5.36633725e-08,  -1.80643006e-07,
  -1.76695163e-07,  -6.02879224e-08,   6.36137567e-08,   1.82948329e-07,
   1.67709606e-07,  -6.95097238e-08,  -2.81543035e-07,  -1.62157223e-07,
   1.59995687e-07,   2.76320615e-07,   7.91191565e-08,  -1.55089035e-07,
  -2.00122815e-07,  -8.37247404e-08,   8.51732622e-08,   2.04441857e-07,
   1.47871244e-07,  -8.92194564e-08,  -2.63290375e-07,  -1.44085068e-07,
   1.43203505e-07,   2.59751664e-07,   9.45666953e-08,  -1.39899410e-07,
  -2.14332489e-07,  -9.76492812e-08,   9.82430663e-08,   2.17206296e-07,
   1.35851448e-07,  -1.00919826e-07,  -2.52237817e-07,  -1.33361007e-07,
   1.33041511e-07,   2.49923380e-07,   1.03909218e-07,  -1.30892598e-07};


extern const double mey512[512]=
{  7.10493973e-01,   4.50105050e-01,  -3.38583052e-03,  -1.49893593e-01,
   3.38163200e-03,   8.97669622e-02,  -3.37437293e-03,  -6.39386112e-02,
   3.36395551e-03,   4.95432853e-02,  -3.35052556e-03,  -4.03450647e-02,
   3.33435448e-03,   3.39457751e-02,  -3.31560344e-03,  -2.92265276e-02,
   3.29420529e-03,   2.55950296e-02,  -3.26998220e-03,  -2.27082463e-02,
   3.24288110e-03,   2.03537102e-02,  -3.21309123e-03,  -1.83929344e-02,
   3.18092593e-03,   1.67319834e-02,  -3.14658664e-03,  -1.53049395e-02,
   3.11004467e-03,   1.40640442e-02,  -3.07115860e-03,  -1.29736823e-02,
   3.02990966e-03,   1.20067463e-02,  -2.98651908e-03,  -1.11424056e-02,
   2.94132983e-03,   1.03645576e-02,  -2.89457042e-03,  -9.66055605e-03,
   2.84623667e-03,   9.02015914e-03,  -2.79620900e-03,  -8.43482325e-03,
   2.74448774e-03,   7.89738929e-03,  -2.69131044e-03,  -7.40199712e-03,
   2.63703363e-03,   6.94396076e-03,  -2.58189657e-03,  -6.51946507e-03,
   2.52590301e-03,   6.12518664e-03,  -2.46893857e-03,  -5.75806376e-03,
   2.41100600e-03,   5.41532586e-03,  -2.35234263e-03,  -5.09465878e-03,
   2.29330211e-03,   4.79426613e-03,  -2.23411832e-03,  -4.51270574e-03,
   2.17478716e-03,   4.24861679e-03,  -2.11518405e-03,  -4.00057150e-03,
   2.05529934e-03,   3.76716783e-03,  -1.99535580e-03,  -3.54724443e-03,
   1.93569056e-03,   3.33998117e-03,  -1.87651914e-03,  -3.14476677e-03,
   1.81781735e-03,   2.96095098e-03,  -1.75943900e-03,  -2.78771644e-03,
   1.70135136e-03,   2.62418785e-03,  -1.64375292e-03,  -2.46966028e-03,
   1.58695542e-03,   2.32371079e-03,  -1.53114803e-03,  -2.18607515e-03,
   1.47627945e-03,   2.05640760e-03,  -1.42217568e-03,  -1.93415902e-03,
   1.36877577e-03,   1.81869136e-03,  -1.31624961e-03,  -1.70951048e-03,
   1.26488015e-03,   1.60638137e-03,  -1.21482769e-03,  -1.50920820e-03,
   1.16601218e-03,   1.41779664e-03,  -1.11823104e-03,  -1.33173436e-03,
   1.07139507e-03,   1.25050747e-03,  -1.02564633e-03,  -1.17373488e-03,
   9.81240476e-04,   1.10128507e-03,  -9.38311185e-04,  -1.03315730e-03,
   8.96752444e-04,   9.69244962e-04,  -8.56336572e-04,  -9.09217036e-04,
   8.16950091e-04,   8.52635196e-04,  -7.78711755e-04,  -7.99188920e-04,
   7.41854891e-04,   7.48812799e-04,  -7.06491882e-04,  -7.01568186e-04,
   6.72496507e-04,   6.57407050e-04,  -6.39621979e-04,  -6.16053707e-04,
   6.07736844e-04,   5.77122285e-04,  -5.76943027e-04,  -5.40352082e-04,
   5.47458177e-04,   5.05725097e-04,  -5.19380158e-04,  -4.73347879e-04,
   4.92569378e-04,   4.43215546e-04,  -4.66766820e-04,  -4.15093655e-04,
   4.41829924e-04,   3.88635799e-04,  -4.17850610e-04,  -3.63619059e-04,
   3.95037601e-04,   3.40061635e-04,  -3.73480883e-04,  -3.18104773e-04,
   3.53034002e-04,   2.97776851e-04,  -3.33432052e-04,  -2.78875312e-04,
   3.14527524e-04,   2.61084311e-04,  -2.96408287e-04,  -2.44210221e-04,
   2.79279866e-04,   2.28299299e-04,  -2.63229856e-04,  -2.13519662e-04,
   2.48110176e-04,   1.99925402e-04,  -2.33655007e-04,  -1.87338559e-04,
   2.19716596e-04,   1.75466793e-04,  -2.06383186e-04,  -1.64138929e-04,
   1.93861253e-04,   1.53422646e-04,  -1.82239870e-04,  -1.43506489e-04,
   1.71372913e-04,   1.34464008e-04,  -1.60996961e-04,  -1.26135762e-04,
   1.50967048e-04,   1.18247020e-04,  -1.41374558e-04,  -1.10643329e-04,
   1.32429417e-04,   1.03408242e-04,  -1.24224420e-04,  -9.67453457e-05,
   1.16617400e-04,   9.07424406e-05,  -1.09349093e-04,  -8.52535646e-05,
   1.02278849e-04,   8.00167265e-05,  -9.55025085e-05,  -7.48895008e-05,
   8.92345538e-05,   6.99667838e-05,  -8.35724159e-05,  -6.54628522e-05,
   7.83786137e-05,   6.14755686e-05,  -7.33985944e-05,  -5.78684340e-05,
   6.84964282e-05,   5.43883476e-05,  -6.37726560e-05,  -5.09012284e-05,
   5.94464334e-05,   4.75097984e-05,  -5.56198142e-05,  -4.44356672e-05,
   5.21598775e-05,   4.17835628e-05,  -4.88165543e-05,  -3.94234083e-05,
   4.54583110e-05,   3.71081051e-05,  -4.21899877e-05,  -3.47091780e-05,
   3.92349341e-05,   3.23345806e-05,  -3.66992861e-05,  -3.02108027e-05,
   3.44540868e-05,   2.84471212e-05,  -3.22531073e-05,  -2.69176971e-05,
   2.99685273e-05,   2.53793769e-05,  -2.77087707e-05,  -2.37073561e-05,
   2.57006403e-05,   2.20130032e-05,  -2.40535937e-05,  -2.05259829e-05,
   2.26418630e-05,   1.93585240e-05,  -2.12222765e-05,  -1.83875304e-05,
   1.96699394e-05,   1.73723962e-05,  -1.80960721e-05,  -1.61906832e-05,
   1.67301470e-05,   1.49559571e-05,  -1.56841670e-05,  -1.38999239e-05,
   1.48347873e-05,   1.31367082e-05,  -1.39411401e-05,  -1.25449747e-05,
   1.28805182e-05,   1.18857527e-05,  -1.17662176e-05,  -1.10381227e-05,
   1.08296779e-05,   1.01170613e-05,  -1.01847652e-05,  -9.35558678e-06,
   9.70989693e-06,   8.86904380e-06,  -9.16586978e-06,  -8.53723207e-06,
   8.43154781e-06,   8.12223683e-06,  -7.62170889e-06,  -7.50412092e-06,
   6.96918995e-06,   6.79877622e-06,  -6.58917319e-06,  -6.24007395e-06,
   6.36131473e-06,   5.94415493e-06,  -6.04757568e-06,  -5.79156150e-06,
   5.52791386e-06,   5.54507187e-06,  -4.91813425e-06,  -5.08539585e-06,
   4.45203824e-06,   4.52903010e-06,  -4.24571307e-06,  -4.11041193e-06,
   4.17966331e-06,   3.94621241e-06,  -4.01664463e-06,  -3.91747061e-06,
   3.63735844e-06,   3.78742988e-06,  -3.15830574e-06,  -3.43723581e-06,
   2.81393946e-06,   2.98379398e-06,  -2.72095664e-06,  -2.66192601e-06,
   2.76043189e-06,   2.58866444e-06,  -2.69565313e-06,  -2.64538810e-06,
   2.40781851e-06,   2.59565951e-06,  -2.01389286e-06,  -2.32092414e-06,
   1.74876250e-06,   1.93836981e-06,  -1.72952932e-06,  -1.68308433e-06,
   1.83764570e-06,   1.67235142e-06,  -1.83675144e-06,  -1.78778666e-06,
   1.60837239e-06,   1.79317530e-06,  -1.26977880e-06,  -1.57017243e-06,
   1.05614207e-06,   1.23616353e-06,  -1.08483019e-06,  -1.02642323e-06,
   1.23754344e-06,   1.05841191e-06,  -1.27815217e-06,  -1.21391182e-06,
   1.08839649e-06,   1.25686512e-06,  -7.85746162e-07,  -1.06907453e-06,
   6.05558756e-07,   7.68064793e-07,  -6.65376007e-07,  -5.89242374e-07,
   8.47059978e-07,   6.50192461e-07,  -9.14631182e-07,  -8.32815083e-07,
   7.49969032e-07,   9.01163121e-07,  -4.70672877e-07,  -7.37143320e-07,
   3.12221293e-07,   4.58378507e-07,  -3.92269027e-07,  -3.00368021e-07,
   5.92783323e-07,   3.80785029e-07,  -6.77882136e-07,  -5.81612514e-07,
   5.29535101e-07,   6.66981139e-07,  -2.65425466e-07,  -5.18870518e-07,
   1.21110269e-07,   2.54972098e-07,  -2.14317680e-07,  -1.10850236e-07,
   4.27083269e-07,   2.04239738e-07,  -5.23588141e-07,  -4.17181631e-07,
   3.85860226e-07,   5.13860803e-07,  -1.31636960e-07,  -3.76307458e-07,
  -3.47381479e-09,   1.22260587e-07,  -9.83037647e-08,   1.26703641e-08,
   3.19049659e-07,   8.92921966e-08,  -4.22985860e-07,  -3.10229518e-07,
   2.92177900e-07,   4.14363873e-07,  -4.43981636e-08,  -2.83760147e-07,
  -8.47125812e-08,   3.61898674e-08,  -2.26525947e-08,   9.27066003e-08,
   2.48602585e-07,   1.48777684e-08,  -3.57385893e-07,  -2.41051851e-07,
   2.31092271e-07,   3.50063432e-07,   1.24833882e-08,  -2.24000884e-07,
  -1.37678980e-07,  -1.93422935e-08,   2.66672827e-08,   1.44304767e-07,
   2.02679507e-07,  -3.30594651e-08,  -3.14627088e-07,  -1.96521276e-07,
   1.91280784e-07,   3.08702428e-07,   4.95500446e-08,  -1.85588018e-07,
  -1.72189419e-07,  -5.50138378e-08,   5.87967660e-08,   1.77427745e-07,
   1.72768179e-07,  -6.38130678e-08,  -2.86782308e-07,  -1.67970564e-07,
   1.65360928e-07,   2.82199598e-07,   7.36773437e-08,  -1.60988361e-07,
  -1.94647503e-07,  -7.78454510e-08,   7.97001447e-08,   1.98617113e-07,
   1.53313254e-07,  -8.34768690e-08,  -2.68676869e-07,  -1.49724187e-07,
   1.48512364e-07,   2.65270067e-07,   8.93557497e-08,  -1.45281733e-07,
  -2.09236571e-07,  -9.24169614e-08,   9.32747632e-08,   2.12135158e-07,
   1.40683752e-07,  -9.60169615e-08,  -2.56927920e-07,  -1.38092286e-07,
   1.37583411e-07,   2.54481563e-07,   9.95215647e-08,  -1.35276034e-07,
  -2.18692260e-07,  -1.01696569e-07,   1.02069320e-07,   2.20741387e-07,
   1.32505104e-07,  -1.03998389e-07,  -2.49323071e-07,  -1.30690960e-07,
   1.30512721e-07,   2.47618852e-07,   1.06095460e-07,  -1.28913016e-07,
  -2.24804194e-07,  -1.07596454e-07,   1.07751334e-07,   2.26212111e-07,
   1.27223603e-07,  -1.09071092e-07,  -2.44414725e-07,  -1.25987801e-07,
   1.25951644e-07,   2.43258848e-07,   1.10333813e-07,  -1.24871296e-07,
  -2.28742719e-07,  -1.11343388e-07,   1.11410960e-07,   2.29686092e-07,
   1.23823810e-07,  -1.12291985e-07,  -2.41257033e-07,  -1.23001982e-07,
   1.23019136e-07,   2.40491414e-07,   1.13057321e-07,  -1.22306368e-07,
  -2.31272322e-07,  -1.13720961e-07,   1.13760301e-07,   2.31890395e-07,
   1.21642470e-07,  -1.14335688e-07,  -2.39232315e-07,  -1.21107556e-07,
   1.21140049e-07,   2.38735804e-07,   1.14801440e-07,  -1.20679489e-07,
  -2.32891432e-07,  -1.15228879e-07,   1.15263272e-07,   2.33288441e-07,
   1.20247808e-07,  -1.15631902e-07,  -2.37938703e-07,  -1.19906138e-07,
   1.19940329e-07,   2.37622682e-07,   1.15914360e-07,  -1.19648238e-07,
  -2.33924171e-07,  -1.16184642e-07,   1.16221580e-07,   2.34174666e-07,
   1.19359023e-07,  -1.16453708e-07,  -2.37114915e-07,  -1.19144431e-07,
   1.19176927e-07,   2.36917106e-07,   1.16622092e-07,  -1.18994732e-07,
  -2.34580653e-07,  -1.16790263e-07,   1.16830525e-07,   2.34736310e-07,
   1.18794590e-07,  -1.16974611e-07,  -2.36592207e-07,  -1.18661688e-07};


// constructors

template<class DataType_t> Meyer<DataType_t>::
Meyer(const Wavelet &w) :
WaveDWT<DataType_t>(w)
{
   pLForward = NULL;
   pLInverse = NULL;
   pHForward = NULL;
   pHInverse = NULL;
   setFilter();
}

template<class DataType_t> Meyer<DataType_t>::
Meyer(const Meyer<DataType_t> &w) :
WaveDWT<DataType_t>(w)
{
   pLForward = NULL;
   pLInverse = NULL;
   pHForward = NULL;
   pHInverse = NULL;
   setFilter();
}

template<class DataType_t> Meyer<DataType_t>::
Meyer(int m, int tree, enum BORDER border) :
WaveDWT<DataType_t>(m,m,tree,border)
{
   pLForward = NULL;
   pLInverse = NULL;
   pHForward = NULL;
   pHInverse = NULL;
   setFilter();
}

// destructor
template<class DataType_t>
Meyer<DataType_t>::~Meyer()
{
   if(pLForward) free(pLForward);
   if(pLInverse) free(pLInverse);
   if(pHForward) free(pHForward);
   if(pHInverse) free(pHInverse);
}

// clone
template<class DataType_t>
Meyer<DataType_t>* Meyer<DataType_t>::Clone() const
{
  return new Meyer<DataType_t>(*this);
}

template<class DataType_t>
int Meyer<DataType_t>::getMaxLevel()
{
   int i;
   int maxLevel = 0;
   int n = int(this->nWWS);
   int m = this->m_L;
   int N[10] = {512,384,256,192,128,96,64,48,32,16};

   if(this->m_TreeType == 2) {                              // adaptive wavelet
      for(i=0; i<10; i++) {if(N[i] == this->m_L/2) break;}  // find initial wavelet

      while((n>=2*m) && !(n&1)){
	 n /= 2;
	 if(i < 10) m = N[i++]*2;
	 else       m = 32;
	 maxLevel++;
      }
      return maxLevel;
   }

   else return Wavelet::getMaxLevel(n);
}

//: calculate wavelet filter
//!param: taper function order n
//!param: beta(n,n) - value of Euler's beta function
//!param: integration step
template<class DataType_t>
double Meyer<DataType_t>::filter(int n, double B, double s)
{
   int i;
   int k = this->m_H;
   double pi3 = PI/3;
   double pi2 = PI/2;
   double x = s;
   double y = 0;
   double z = 0;
   double v = 0;

   if(n<2) n = 2;
   B = pow(B,1./(n-1));

   double* H = (double*)malloc(k*sizeof(double));
   for(i=0; i<k; i++) H[i] = 0.;

   while(x<1.){
      v += s*pow(x*(1-x)/B,n-1);
      y  = cos(v*pi2);
      z  = pi3*(1+x);
      for(i=0; i<k; i++) H[i] += y*cos(i*z)*s;
      x += s;
//      if(int(x/0.1)*0.1==x) printf(" %6.5f \n",x);
   }

   v = y = 0.;
   printf("\n");
   for(i=0; i<k; i++){
      H[i] += i ? sin(pi3*i)/pi3/i : 1;
      H[i] *= sqrt(2.)/3.;
      v += 2*H[i]*H[i];
      y += 2*H[i];
      printf(" %16.8e,",H[i]);
      if(((i+1)/4)*4==i+1) printf("\n");
   }
   printf("\n normalization: %9.6f, orthogonality: %9.6f \n",y-H[0],v-H[0]*H[0]);
   free(H);
   return v-H[0]*H[0];
}


template<class DataType_t>
void Meyer<DataType_t>::setFilter()
{
   const double* pF;
   this->m_H = (this->m_H>>1)<<1;
   int n = this->m_H/2;
   switch(n)
   {
      case 512: pF =  mey512; break;
      case 384: pF =  mey384; break;
      case 256: pF =  mey256; break;
      case 192: pF =  mey128; break;
      case 128: pF =  mey064; break;
      case  96: pF =  mey032; break;
      case  64: pF =  mey016; break;
      case  48: pF =  mey008; break;
      case  32: pF =  mey004; break;
      case  16: pF =  mey;    break;
      default:  pF =  mey;
         this->m_H = this->m_L = 32;
         break;
   }

   if(!pLForward) {
      pLInverse = (double*)malloc(this->m_H*sizeof(double));
      pLForward = (double*)malloc(this->m_H*sizeof(double));
      pHInverse = (double*)malloc(this->m_H*sizeof(double));
      pHForward = (double*)malloc(this->m_H*sizeof(double));
   }
   else {
      pLInverse = (double*)realloc(pLInverse,this->m_H*sizeof(double));
      pLForward = (double*)realloc(pLForward,this->m_H*sizeof(double));
      pHInverse = (double*)realloc(pHInverse,this->m_H*sizeof(double));
      pHForward = (double*)realloc(pHForward,this->m_H*sizeof(double));
   }

   double* temp = (double*)malloc(this->m_H*sizeof(double));

   temp[0] = 0.;
   for(int i=0; i<this->m_H/2; i++){
      temp[i+this->m_H/2] = pF[i];
      temp[this->m_H/2-i] = pF[i];
   }

//  LP filter for db3:  h0  h1  h2  h3  h4  h5
//  HP filter for db3:  h5 -h4  h3 -h2  h1 -h0
// iLP filter for db3:  h4  h1  h2  h3  h0  h5
// iHP filter for db3:  h5 -h0  h3 -h2  h1 -h4

   for(int i=0; i<this->m_H; i+=2){

      pLForward[i]   = temp[i];
      pLForward[i+1] = temp[i+1];
      pHForward[i]   = temp[this->m_H-1-i];
      pHForward[i+1] = -temp[this->m_H-2-i];

      pLInverse[i]   = temp[this->m_H-2-i];
      pLInverse[i+1] = temp[i+1];
      pHInverse[i]   = temp[this->m_H-1-i];
      pHInverse[i+1] = -temp[i];

   }

   this->m_WaveType = MEYER;
   free(temp);
}

// forward function does one step of forward transformation.
// <level> input parameter is the level to be reconstructed
// <layer> input parameter is the layer to be reconstructed.
template<class DataType_t>
void Meyer<DataType_t>::forward(int level,int layer)
{
   if(this->m_TreeType == 2 && !layer){                              // adaptive wavelet
      int i;
      int N[10] = {512,384,256,192,128,96,64,48,32,16};
      for(i=0; i<10; i++) {if(N[i] == this->m_L/2) break;}  // find initial wavelet
      if(i+level > 9) { this->m_H = 32; }
      else { this->m_H = N[i+level]*2; }
      setFilter();
   }
   this->forwardFWT(level, layer, pLForward, pHForward);
}

// inverse function does one step of inverse transformation.
// <level> input parameter is the level to be reconstructed
// <layer> input parameter is the layer to be reconstructed.
template<class DataType_t>
void Meyer<DataType_t>::inverse(int level,int layer)
{
   if(this->m_TreeType == 2 && !layer) {                   // adaptive wavelet
      int i;
      int N[10] = {512,384,256,192,128,96,64,48,32,16};
      for(i=0; i<10; i++) {if(N[i] == this->m_L/2) break;}  // find initial wavelet
      if(i+level > 9) { this->m_H = 32; }
      else { this->m_H = N[i+level]*2; }
      setFilter();
   }
   this->inverseFWT(level, layer, pLInverse, pHInverse);
}

// instantiations

#define CLASS_INSTANTIATION(class_) template class Meyer< class_ >;

CLASS_INSTANTIATION(float)
CLASS_INSTANTIATION(double)
//CLASS_INSTANTIATION(std::complex<float>)
//CLASS_INSTANTIATION(std::complex<double>)

#undef CLASS_INSTANTIATION

//}  // end namespace wat
//}  // end namespace datacondAPI

