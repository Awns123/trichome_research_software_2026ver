# 알려진 한계

## 1. 작은 분류군 표본과 행렬 자료

주 분석의 추론 단위는 여섯 taxon label입니다. 거리행렬 상삼각의 15개 원소는 서로 독립된 15개 표본이 아닙니다. `6!=720`개 taxon-label 전수순열은 이 구조에 맞는 선택이지만, 여섯 종만으로 효과의 안정성과 일반성을 평가하기에는 검정력이 작습니다.

*B. suaveolens* 제외 시 원 Spearman 상관이 `0.775`에서 `0.333`으로 감소하고 유의하지 않았습니다. 이는 종 구성 민감성이지 특정 deep split의 증명이나 외부검증이 아닙니다.

## 2. 형질 블록의 탐색적 선택

주 형질은 `tortuosity`, `curvature_mean_rad`, `log1p(curvature_sum_rad / length_um)`입니다. 여러 단일형질과 블록을 검토한 뒤 선택됐으며 사전등록된 주 분석이 아닙니다. 따라서 원 `p=0.044`는 탐색적 결과로 읽어야 합니다.

본문의 Curvature per Length 설명과 주코드의 `log1p` 변환 정의도 완전히 일치하지 않습니다. 저장소에서는 실제 주코드 정의를 우선해 명시합니다.

## 3. 반복단위와 기반자료

분석표는 165개의 image-level median으로 구성됩니다.

| 종 | image median 수 |
|---|---:|
| *B. suaveolens* | 30 |
| *C. annuum* var. *grossum* | 31 |
| *L. chinense* | 31 |
| *P. × hybrida* ‘Dreams Red’ | 18 |
| *P. alkekengi* var. *franchetii* | 28 |
| *S. nigrum* | 27 |

plant·leaf·batch 식별자가 없어 독립 식물 개체 수준의 생물학적 반복과 촬영·잎 내부 유사반복을 분리할 수 없습니다. 자료 생성 일부는 2025년 팀 연구에서 이루어졌습니다.

## 4. 형태 측정의 민감성

곡률 형질은 raster skeleton, 화소 계단 현상, 중심선 방향과 회전에 민감합니다. 2026년 후속 회전 감사에서는 `tortuosity`, mean curvature, curvature/length 중앙값 변화가 각각 약 4.13%, 23.98%, 30.77%였습니다. 이는 HOBY 당시 결과가 아니라 이후 감사 결과입니다.

원 image-level CSV에는 비초점 필드의 비정상값도 있습니다. 주 분석이 사용한 세 형질 외의 `contraction_ratio`, perimeter, ellipse 계열 값을 검증된 측정치로 재사용하면 안 됩니다.

## 5. 분자입력과 tree 재현성

역사적 아카이브에는 원 FASTA에서 보관 tree까지 정확히 이어지는 MSA, supermatrix, partition, 명령, 프로그램 버전, model, seed와 실행 로그가 완전하게 남아 있지 않습니다. 따라서 역사적 공개 코드는 **처리된 분자거리 행렬 이후의 Mantel 계산**은 재현하지만 원 tree 추론 전체를 재현하지 않습니다.

사후 감사에서 *Brugmansia* strand 방향 오류와 전체 자료의 종 식별·partial-window·raw provenance 문제가 확인됐습니다. 역사적 분자거리 행렬은 원 결과 재현용으로만 보존하며 현재의 검증된 입력으로 사용하지 않습니다.

## 6. 통계 코드 한계

- `solanaceae_permanova_permdisp_code.py`의 PERMDISP 순열은 label을 바꿀 때 group centroid를 새로 계산하지 않습니다. 해당 `p=0.045`는 격리된 역사적 출력이며 신뢰 결과로 사용하지 않습니다.
- Ward bootstrap support는 낮고, 형질 열 재표집은 상관·파생 형질의 비독립성을 충분히 반영하지 않습니다.
- exact Mantel 자체가 입력 행렬의 오류나 사후 형질 선택을 교정해 주지는 않습니다.
- 일부 역사적 코드에는 random seed와 package version이 고정돼 있지 않습니다.

## 7. 코드·도표 재현성

- 원 분석·도표 코드 다수는 `/mnt/data`를 고정 경로로 사용합니다.
- 원 notebook 부재로 여덟 도표 코드는 저장 산출물에서 사후 재구성됐습니다.
- 일부 그림은 저장 수치를 직접 사용하며, Figure 6 코드는 재생성 CSV 파일명이 맞지 않으면 manuscript fallback 값을 사용합니다.
- `requirements.txt`는 필요한 package 이름만 제시합니다. 2026년 당시의 정확한 환경을 소급해 고정하지 않습니다.

2026-05-05 휴대용 재계산 코드를 현재 공개자료에 실행하면 Spearman 주결과는 `r_M=0.775`, `p=0.044444`로 회수되지만 Pearson은 보관 요약의 `0.865307` 대신 `0.865481`을 냅니다. 재계산 형태거리 행렬과 보관 행렬의 최대 절대차는 약 `0.001862`입니다. 또한 contraction 민감도 Pearson p-value는 보관 요약 `0.5722`와 달리 `0.197222`입니다. 따라서 휴대용 코드는 원 분석의 완전한 canonical 재현본이 아니라 사후 재구성으로 취급합니다.

## 8. 허용되지 않는 주장

이 자료로 “세 형질만으로 가지과 종을 안정적으로 분류했다”, “모용 형태로 계통을 복원했다”, “일반적인 형태–계통 대응을 입증했다” 또는 “교정 분석이 계통 신호의 부재를 증명했다”고 주장하지 않습니다.
