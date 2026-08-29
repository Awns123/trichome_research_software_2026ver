# 출처계보와 파일 경계

## 원본 수집 스냅샷

- 논리적 파일명: `2026.zip`
- 크기: `815,716,192 bytes`
- SHA-256: `d35caf7593a3a343c5f6f7c234490ce761a4a44517a0f95366e5296fc7eb53cd`
- ZIP entry 수: `505`
- 해시 재확인: 2026-08-26

외장하드의 실제 경로는 공개하지 않습니다. 이 해시는 현재 보관 중인 수집 ZIP의 바이트 식별자일 뿐, 각 파일의 저자성·작성시점·과학적 타당성 또는 아카이브의 완전성을 독립적으로 증명하지 않습니다.

ZIP 내부 수정시각도 출처 보조정보입니다. Git 개발이력이나 제3자 timestamp 인증이 아니며, 공개 Git commit을 2026년 연구 시점으로 소급하지 않습니다.

## 공개 층위

### 1. 역사적 분석

`historical_2026/04_code/morphology_phylogeny_analysis/`의 네 `solanaceae_*.py` 파일은 2026-04-12 HOBY 분석 세션에 가장 가까운 보존 코드입니다. 원본의 `/mnt/data` 고정 경로와 구현상 한계를 숨기지 않기 위해 내용을 수정하지 않았습니다.

같은 폴더의 `common_utils.py`와 `figure_generation/`의 여덟 파일은 원 notebook이 남아 있지 않은 상태에서 저장된 CSV·그림을 바탕으로 2026년 4월 재구성한 도표 코드입니다. 원 notebook의 역사적 원본이라고 주장하지 않습니다.

`portable_reproduction/`의 두 파일은 2026-05-05에 정리된 휴대용 사후 재계산 코드입니다. 역사적 처리자료에서 Mantel 수치와 일부 그림을 재생성하지만, raw FASTA부터 기존 tree를 다시 만드는 pipeline은 아닙니다.

### 2. 제출 뒤 독립 감사

`post_submission_audit/`는 원 ZIP 밖에서 2026년 7월 진행된 Codex 보조 독립 재분석에서 선별한 소규모 자료입니다. 원 연구 당시 학생 작업으로 소급하지 않습니다. 저장 snapshot에는 학생의 독립 전체 재실행이 완료되지 않았다고 기록돼 있습니다.

공개본에는 focal 여섯 종의 accession 감사표, 선별 sequence set, alignment, tree와 거리행렬만 포함했습니다. 사용자 로컬 경로가 든 실행 JSON, 다운로드 cache, MAFFT·IQ-TREE 바이너리, 전체 NCBI record와 17종 원자료는 제외했습니다.

### 3. 2026-08-26 공개 문서와 검산기

루트의 설명 문서, `.gitignore`, `requirements.txt`, `SHA256SUMS.txt`, `PROVENANCE_FILES.csv`와 `post_submission_audit/code/verify_mantel.py`는 공개 준비 과정에서 AI 보조로 작성·정리했습니다. 역사적 학생 산출물로 제시하지 않습니다.

## 제외한 자료

- 원본·중첩 ZIP과 전체 파일목록
- 논문 PDF, 서명된 AI 공개서, 상장, 심사평
- 연구일지, HWPX, DOCX, XLSX와 비공개 과정증거
- 원본 SEM 이미지와 팀 원자료
- 2025 코드와 가상환경
- 개인정보 또는 로컬 사용자 경로가 포함된 실행 JSON·로그
- MAFFT, IQ-TREE, BLAST 등의 바이너리와 vendored dependency
- SlHD8·orthology 후속 감사 코드

각 공개 파일의 SHA-256은 `SHA256SUMS.txt`에서, 원 ZIP 내부 수정시각과 층위는 `PROVENANCE_FILES.csv`에서 확인할 수 있습니다.

