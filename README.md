# 2026 모용 형태–분자계통 거리 연구 기록

이 저장소는 2026년 학생 연구 **“Correspondence Between Trichome Morphology and Phylogenetic Distance in Six Solanaceae Species”**의 역사적 분석 코드·처리자료와, 제출 뒤 이루어진 분자입력 사후감사를 서로 구분해 보존합니다.

> 연구 시점: 2026년  
> 공개 준비일: 2026-08-26  
> 저장소 성격: 역사적 분석 스냅샷과 명시적으로 분리된 사후 교정 기록  
> 현재 과학적 상태: 원래 보고한 양의 형태–계통 대응 결론은 유지하지 않음

먼저 [CORRECTION_NOTICE.md](CORRECTION_NOTICE.md)를 읽어 주세요. 이 저장소의 역사적 입력·결과는 오류가 교정된 현재 결론이 아니라, 무엇이 바뀌었는지를 추적하기 위한 보존 자료입니다.

## 1분 요약

| 질문 | 현재 답 |
|---|---|
| 무엇을 물었나? | 가지과 여섯 종에서 세 모용 형질로 계산한 형태거리와 `rbcL+matK` 계통거리 사이에 대응이 있는지 물었습니다. |
| 당시 결과는? | 보관된 분자거리 행렬을 사용한 Spearman exact Mantel은 `r_M=0.775`, 단측 `p=0.044`; Pearson은 `r_M=0.865`, `p=0.046`이었습니다. 여섯 종 라벨의 `6!=720`개 순열을 모두 사용했습니다. |
| 당시에도 불안정했나? | *B. suaveolens*를 제외하면 Spearman `r_M=0.333`, 논문에 보고된 양측 `p=0.575`로 신호가 사라졌습니다. |
| 뒤에 무엇이 확인됐나? | 제출 뒤 독립 감사에서 보관된 *B. suaveolens* `matK` 입력의 strand 방향 문제가 확인되었습니다. 기존 accession을 주석된 5′→3′ 방향으로 다시 추출한 분석은 Spearman `r_M=-0.243`, 단측 `p=0.676`이었습니다. |
| 현재 결론은? | 원래의 양의 대응을 현재의 생물학적 결론으로 유지하지 않습니다. 다만 여섯 종·두 plastid marker 결과만으로 모용에 계통 신호가 없다고 결론 내릴 수도 없습니다. |
| 이 저장소의 가치는? | 유의한 수치를 지키는 대신, 원 분석·사후 복원·독립 감사를 구분하고 입력 오류가 결론을 바꾸는 과정을 추적 가능하게 만든 데 있습니다. |

정확한 결과 층위와 허용되는 해석은 [SCIENTIFIC_STATUS.md](SCIENTIFIC_STATUS.md)에 정리했습니다.

## 저장소 구조

```text
.
├─ historical_2026/
│  ├─ 03_processed_data/       # 당시 분석에 사용된 소규모 처리자료
│  ├─ 04_code/
│  │  ├─ morphology_phylogeny_analysis/  # 당시 분석 세션에 가까운 코드
│  │  ├─ figure_generation/             # 저장 산출물에서 사후 재구성한 도표 코드
│  │  └─ portable_reproduction/         # 2026-05-05 휴대용 재계산 코드
│  └─ 05_results/              # 보관된 역사적 표·그림·요약
├─ post_submission_audit/
│  ├─ data/                    # 원 행렬, 방향 교정 행렬, accession 감사표
│  ├─ sequences/               # 여섯 종 focal sequence set
│  ├─ trees/                   # 선별된 alignment와 tree 산출물
│  ├─ figures/                 # 감사 결과 비교 그림
│  └─ code/                    # 공개 준비용 결정론적 수치 검산기
├─ PROVENANCE.md
├─ AI_AND_CONTRIBUTIONS.md
├─ SCIENTIFIC_STATUS.md
├─ KNOWN_LIMITATIONS.md
├─ VALIDATION.md
└─ SHA256SUMS.txt
```

원본 SEM 이미지, 연구일지, 서명 문서, 논문 PDF, 팀 자료, 가상환경과 실행 바이너리는 개인정보·공동연구 권리·용량 문제 때문에 포함하지 않았습니다. 2025년 영상 정량화 프로그램도 2025 코드 저장소와 중복되므로 이 저장소에서 제외했습니다.

## 빠른 검산

Python 3 가상환경을 만든 뒤 필요한 패키지를 설치합니다.

```bash
python -m venv .venv

# Windows PowerShell
.venv\Scripts\python -m pip install -r requirements.txt

# macOS/Linux
.venv/bin/python -m pip install -r requirements.txt
```

역사적 처리자료에서 2026-05-05 휴대용 재구성 경로와 그림을 다시 실행합니다.

```bash
python historical_2026/04_code/portable_reproduction/reproduce_main_analysis.py
python historical_2026/04_code/portable_reproduction/reproduce_figures.py
```

원 행렬과 교정 행렬의 exact Mantel 결과를 다시 계산하고 보관된 감사표와 대조합니다.

```bash
python post_submission_audit/code/verify_mantel.py --check
```

첫 두 명령은 **보관된 잘못된 분자입력을 사용하는 사후 재구성 경로**입니다. 원 Spearman 수치는 회수하지만 일부 Pearson·민감도 수치는 보관 요약과 일치하지 않습니다. 정확한 차이는 [VALIDATION.md](VALIDATION.md)에 기록했습니다. 세 번째 명령은 보관 행렬 자체에서 원 결과와 교정 결과를 함께 검산합니다.

## 저자성·AI 사용 경계

- 2026 최종 논문은 단독 저자로 제출됐지만, 기반 SEM 자료는 2025년 팀 연구에서 생성됐습니다.
- 서명된 AI 사용 공개서에 따르면 연구 질문·분석 설계·해석과 초기 코드 설계는 연구자가 수행했고, 생성형 AI는 일부 세부 구현과 영문 번역·편집을 보조했으며, 연구자가 최종 코드를 검토·검증했습니다.
- 2026년 7월의 분자입력 감사는 Codex 보조 독립 재분석입니다. 보존된 스냅샷 시점에는 학생이 전체 파이프라인을 독립 재실행하지 않았습니다.
- 이 GitHub 저장소와 commit 시각은 과거의 개인 저자성이나 기여율을 독립적으로 증명하지 않습니다.

자세한 내용은 [AI_AND_CONTRIBUTIONS.md](AI_AND_CONTRIBUTIONS.md)와 [PROVENANCE.md](PROVENANCE.md)를 확인해 주세요.

## 이용 조건

현재 별도의 오픈소스 라이선스를 부여하지 않았습니다. 공개 열람은 가능하지만 명시적인 허가 없이 복제·수정·재배포할 권리를 부여한다는 뜻은 아닙니다. 공동연구 자료와 외부 데이터의 권리를 확인한 뒤 라이선스를 별도로 정할 수 있습니다.
