0. 코드가 저장되어 있는 폴더(기준 폴더) 내에 있는 fits파일을 사용하여 코드가 돌아가고,
   기준 폴더 내에 'results' 폴더와 'fits파일 이름'의 폴더가 생성되고 그 안에 결과들이 저장됩니다.
   00_~.py 부터 순서대로 각 단계에 맞는 결과를 얻고, 그 결과를 이용하여 다음 코드를 돌리는 방식입니다.
   함수는 일단 따로 빼지 않고 만들었습니다. 테스트 중 이었어서 parameter가 정리되지 않고 지저분한 감이 있어요.^^;
   결과 값을 구한 다음 그래프로 출력은 간단히 따로 했었어서 코드에 그 부분이 없을 수 있어요. 


1. stddev vs. width
 -stddev (standard deviation)
   : fitting 결과로 얻은 model과 fitting에 사용된 data point들의 차이 (즉, model과 data의 차이를 나타냅니다.). 일반적인 표준편차.

 -width
   : Gaussian 함수의 FWHM을 나타냅니다. 일반적으로 Gaussian 함수에서 사용되는 sigma와 다음과 같은 관계가 있습니다.
   width(FWHM) = sigma/(2*sqrt(ln(2)))

  *문서와 코드내에서 헷갈리지 않게 구분하려고 했으나, 수정이 안 된 부분이 많이 있을거에요.^^;


2. OLD 버전, NEW 버전
   OLD 버전은 이미지를 회전시키지 않고, 단면의 data point를 얻어서 피팅 한 것이고,
  NEW 버전은 먼저 이미지를 회전시킨 후 피팅 한 것입니다.
   두 버전은 아래 두 개의 코드를 공유합니다.
   	00_m0_center_2018xxxx.py
	01_m1_tail_center_2018xxxx.py
  NEW 버전을 만들면서 00과 01에 조금 손을 대어서 OLD 버전을 돌렸을 때 오류가 날지도 몰라요.^^;

   정리하면, OLD 버전은
	00_m0_center_2018xxxx.py
	01_m1_tail_center_2018xxxx.py
	02_DL_line_profiles_18xxxx.py
	03_spectrum_18xxxx.py   
  이고,
   NEW 버전은
	00_m0_center_2018xxxx.py
	01_m1_tail_center_2018xxxx.py
	02_image_rotation_2018xxxx.py
	03_fitting_2018xxxx.py
  입니다.


3. image rotation에서 x또는 y 축으로의 평행이동.

  회전 전의 이미지의 각 픽셀들이 제 1 사분면에 있다고 할 때 (x, y좌표가 모두 양수),
 x축 기준 +90도 내로 회전할 때 일부 픽셀들의 x 좌표가 음수가 되고,
 x축 기준 -90도 내로 회전할 때 일부 픽셀들의 y 좌표가 음수가 됩니다.
  회전한 이미지array에서 음수가 된 index를 양수가 되도록 적당히 평행이동시켜서 사용했는데,
 구해진 각도에 따라서 갈래를 나누어 처리하는 부분이 필요할 것입니다.
  
 
4. fits파일 (data 폴더)
  수현이가 주었던 H alpha 라인 있는 타겟의 fits파일이 있습니다. 
 모두 100-c이고, 리덕션이 되어있어요. BS_로 시작하는 파일을 사용하면 됩니다.
  EOS에서 다운받은 GD71 스펙트럼 data도 있어요.

5. 참고자료 (ref 폴더)
  오성아씨와 수현이의 발표 자료와 논문이 있습니다. 논문은 읽어보지는 않았는데 참고하시면 좋을 것 같아요.
