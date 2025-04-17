# ProteinNet

ProteinNet은 단백질 구조에 대한 머신러닝을 위한 표준화된 데이터셋입니다.
단백질 시퀀스와 구조, 다중서열정렬(MSA), PSSMs, 그리고 표준화된 훈련용, 평가용, 테스트용 데이터를 제공합니다.
ProteinNet은 단백질 구조 예측 대회(CASP)를 기반으로 만들어 졌습니다.
CASP는 이미 알고는 있지만 공개되지 않은, 단백질 구조에 대한 블라인드 테스트를 수행하는 대회입니다.
데이터가 많거나 적은 환경안에서도 새로운 방법을 평가하기 위해서, 다양한 데이터 셋의 크기를 하나로 이어지는 데이터 형태로 제공합니다.

** 이건 아직 개발중에 있습니다 **

데이터 셋으로 만들기 위한 원본 데이터들은 아직 이용가능하지 않습니다. 하지만 ProteinNet 12에 쓰인 raw MSA data (4TB)는 요청에 의해 제공 될 수 있습니다.
더 자세한 내용은 [여기](https://github.com/aqlaboratory/proteinnet/blob/master/docs/raw_data.md)를 클릭하세요

### 동기

단백질 구조 예측은 생화학 분야에서 가장 어려운 문제 중 하나 입니다. 이 문제는 생물학과 화학분야에서 중요한 주제이지만 머신러닝 커뮤니티에서는 생소한 분야입니다.
이는 두가지 이유때문이라고 추측됩니다. 1. 높은 진입 장벽 2. 표준화의 부재 이 두가지 문제가 해결된다면 단백질 구조 예측은 비전인식, 음성인식과 더불어 머신러닝의 주요 분야가 될 수 있습니다.
[ImageNet](http://www.image-net.org)이 컴퓨터 비전 기술 발전의 원동력이 되었듯이 ProteinNet은 머신 러닝 분야의 단백질 구조 부분에서 누구든 쉽게 시작 할 수 있도록 표준화된 데이터셋과 트레이닝, 평가, 테스트를 제공 할 것입니다.


### 접근법

CASP 대회는 2년에 한번 열립니다. 이 대회에서는 최근에 밣혀 졌지만, 아직 공개되지 않은 단백질 서열에 대한 구조를 전세계 참가자들이 해결하게 됩니다.
대회 참가자들은 이런 구조들에 대해 블라인드 예측을 하고 정확성을 평가받게 됍니다. 따라서 CASP 구조는 특정 시점에서 얼마나 예측이 잘 되었는가에 대한 표준화된 기준점을 제공합니다. ProteinNet의 기본적인 생각은 CASP 테스트 셋을 사용하여 CASP에 편승하는 것 입니다. Proteinnet은 훈련, 평가용 데이터를 CASP 실험 이전의 조건을 재설정 함으로써 테스트 셋을 보완합니다. 특히 Proteinnet은 사용 가능한 서열과 구조를 시작 전에 제한합니다. 이건 [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi)와 같은 표준 데이터베이스는 히스토리 버전을 유지하지 않으므로 중요합니다.
우리는 [UniParc](http://www.uniprot.org/uniparc/)의 타임리셋 버전과 [JGI](https://img.jgi.doe.gov/)에서 metagenomic 시퀀스를 사용하는데, MSA를 도출하는 시퀀스 데이터베이스 구축을 위해 이 두가지를 사용합니다.
더 나아가 Proteinnet은 쉬운것부터 어려운 것까지 세분화된 평가 데이터를 제공합니다.
쉬운 난이도에서는 모델이 단백질 구조의 마이너한 변화(이를테면 돌연변이)를 예측하는 능력이 어느정도 되는지 평가하는데 유용합니다.
어려운 난이도에서는 모델이 완전히 새로운 단백질 접힘(CASP Free Modeling)을 예측하는데 도움이 됩니다.
이런 평가 데이터는 모델이 데이터셋의 분포 변화를 얼마나 잘 커버하는지 테스트하기 위한 가반성 문제를 제공합니다.
우리는 이런 점을 Proteinnet의 가장 어려운 평가 셋이 CASP FM보다 어렵다는 점에서 알 수 있었습니다.

### 다운로드

Proteinnet의 기록은 두가지 형태로 제공됩니다. 하나는 사람과 기계 모두 읽을 수 있는 텍스트 파일(프로그래밍 가능한 파일), 다른 하나는 텐서플로에 특화된 TFRecord파일입니다. 파일 형식에 대한 더 많은 정보는 [여기](https://github.com/aqlaboratory/proteinnet/blob/master/docs/proteinnet_records.md#file-formats)를 클릭하세요.

| CASP7 | CASP8 | CASP9 | CASP10 | CASP11 | CASP12* |
| --- | --- | --- | --- | --- | --- |
| [Text-based](https://sharehost.hms.harvard.edu/sysbio/alquraishi/proteinnet/human_readable/casp7.tar.gz) | [Text-based](https://sharehost.hms.harvard.edu/sysbio/alquraishi/proteinnet/human_readable/casp8.tar.gz) | [Text-based](https://sharehost.hms.harvard.edu/sysbio/alquraishi/proteinnet/human_readable/casp9.tar.gz) | [Text-based](https://sharehost.hms.harvard.edu/sysbio/alquraishi/proteinnet/human_readable/casp10.tar.gz) | [Text-based](https://sharehost.hms.harvard.edu/sysbio/alquraishi/proteinnet/human_readable/casp11.tar.gz) | [Text-based](https://sharehost.hms.harvard.edu/sysbio/alquraishi/proteinnet/human_readable/casp12.tar.gz) |
| [TF Records](https://sharehost.hms.harvard.edu/sysbio/alquraishi/proteinnet/tfrecords/casp7.tar.gz) | [TF Records](https://sharehost.hms.harvard.edu/sysbio/alquraishi/proteinnet/tfrecords/casp8.tar.gz) | [TF Records](https://sharehost.hms.harvard.edu/sysbio/alquraishi/proteinnet/tfrecords/casp9.tar.gz) | [TF Records](https://sharehost.hms.harvard.edu/sysbio/alquraishi/proteinnet/tfrecords/casp10.tar.gz) | [TF Records](https://sharehost.hms.harvard.edu/sysbio/alquraishi/proteinnet/tfrecords/casp11.tar.gz) | [TF Records](https://sharehost.hms.harvard.edu/sysbio/alquraishi/proteinnet/tfrecords/casp12.tar.gz) |

* CASP 12 테스트 셋은 미완성입니다.(엠바고중임) 엠바고 끝나면 공개하겠습니다.

### 문서
* [ProteinNet Records](docs/proteinnet_records.md)
* [Splitting Methodology](docs/splitting_methodology.md)
* [Raw Data](docs/raw_data.md)
* [FAQ](docs/FAQ.md)

### PyTorch Parser
Proteinnet은 텐서플로기반 공식 파서를 제공합니다. 파이토치기반의 파서는 [Jeppe Hallgren](https://github.com/JeppeHallgren)씨가 만들었고, [여기](https://github.com/OpenProtein/openprotein/blob/master/preprocessing.py)서 이용 할 수 있습니다.

### 인용
인용은 [여기](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2932-0)서 부탁드립니다.(BMC Bioinformatics 링크입니다)

### 감사의 말
이렇게 데이터 셋을 만드 수 있었던 것은 전부 [HMS Laboratory of Systems Pharmacology](http://hits.harvard.edu/the-program/laboratory-of-systems-pharmacology/about/), the [Harvard Program in Therapeutic Science](http://hits.harvard.edu/the-program/program-in-regulatory-science/about/), 그리고 [Harvard Medical School](https://hms.harvard.edu)의 the [Research Computing](https://rc.hms.harvard.edu) 그룹 덕 입니다. 그리고 [Martin Steinegger](https://github.com/martin-steinegger)와 [Milot Mirdita](https://github.com/milot-mirdita)에게도 MMseqs2, HHblits software packages에 대한 많은 도움에 역시 감사를 표합니다. [Sergey Ovchinnikov](http://site.solab.org/)에게는 metagenomic sequences 제공에 대한 감사를 표합니다. [Andriy Kryshtafovych](http://predictioncenter.org/people/kryshtafovych/index.cgi)에게는 CASP 데이터에 대한 도움에 감사를 표합니다. 또 [Sean Eddy](https://github.com/cryptogenomicon)에게는 HMMer software package에 대한 도움을 받아 이에 감사를 표합니다.
이 데이터 셋은 전부 하버드 대학의 the [HMS Research Information Technology Solutions](https://rits.hms.harvard.edu) 그룹이 주도했습니다.

### 펀딩
이 프로젝트는 NIGMS grant P50GM107618 and NCI grant U54-CA225088

# translate to korean by Bue-Von-hon
hoping to be helpful....
