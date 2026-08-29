# 공개본 검증 기록

검증일: 2026-08-26  
운영체제: Windows  
Python: 3.12.13  
상태: 공개 후보 검산 완료

## 검증 환경

```text
numpy==2.5.2
pandas==3.0.5
scipy==1.18.1
scikit-learn==1.9.0
matplotlib==3.11.1
biopython==1.88
```

이 버전은 2026-08-26 공개 검산 환경이며 역사적 환경이 아닙니다. `requirements.txt`에는 당시 정확한 version 기록이 없으므로 package 이름만 적었습니다.

## 통과한 검사

| 검사 | 결과 |
|---|---|
| 공개 Python 구문 분석 | 16/16 통과 |
| `historical_2026` 선별 파일과 원 재현 패키지 SHA-256 대조 | 24/24 일치 |
| `post_submission_audit` 선별 파일과 2026-07-31 감사 snapshot SHA-256 대조 | 44/44 일치 |
| 휴대용 주 분석 실행 | 종료코드 0, CSV 4개 생성 |
| 휴대용 그림 실행 | 종료코드 0, 비어 있지 않은 PNG 3개 생성 |
| 사후감사 exact Mantel 검산 | `CHECK_OK`, 다섯 행렬의 Spearman·Pearson 단측/양측 값 모두 저장표와 일치 |
| focal image-level 자료 수 | 165행; 종별 30/31/31/18/28/27 확인 |

## 휴대용 재구성에서 확인된 불일치

`historical_2026/04_code/portable_reproduction/reproduce_main_analysis.py`는 원 Spearman 주결과 `r_M=0.775000`, `p=0.044444`를 회수했습니다. 그러나 다음 값은 보관 요약과 달랐습니다.

| 항목 | 보관 요약 | 휴대용 재실행 |
|---|---:|---:|
| 주결과 Pearson `r_M` | `0.865307` | `0.865481` |
| focal+contraction Pearson `p` | `0.5722` | `0.197222` |

재계산 형태거리 행렬과 보관 형태거리 행렬의 최대 절대차는 약 `0.001862`였습니다. 따라서 휴대용 코드를 원 분석의 완전한 canonical pipeline으로 판정하지 않았습니다. 보관 행렬 자체를 사용하는 `post_submission_audit/code/verify_mantel.py --check`는 원·교정 수치를 모두 정확히 재계산했습니다.

## 실행 중 경고

Matplotlib은 사용자 cache 폴더에 font cache lock을 저장하지 못했다는 권한 경고를 냈지만, 세 PNG 저장은 완료됐습니다. 저장소 내부 `.mplconfig`를 지정하면 이 환경 의존 경고를 피할 수 있습니다.

## 공개 전 정적 검사

개인 사용자 이름, 학교명, 이메일, 토큰·비밀번호·API key 패턴을 공개 후보에서 검색합니다. 의도적으로 공개한 논문 저자명 외의 개인정보와 사용자 전용 `C:\Users\<name>` 경로는 포함하지 않습니다. IQ-TREE 요약의 `C:\Users\Public` 임시경로와 역사적 코드의 `/mnt/data`는 개인정보가 아니며 출처·한계로 문서화했습니다.

이 검증은 공개본의 파일 무결성과 계산 경로를 확인합니다. 학생의 역사적 저자성, 원 tree의 과학적 타당성, 2026년 당시 실행환경 또는 학생의 독립적인 7월 감사 재실행을 인증하지 않습니다.
