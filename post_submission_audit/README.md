# Post-submission independent audit

이 폴더는 2026년 7월의 Codex 보조 독립 재분석에서 선별한 자료입니다. HOBY 원 연구 당시 학생 개발 코드나 학생 단독 재실행 결과로 제시하지 않습니다. 보존 snapshot 시점에는 학생의 전체 pipeline 독립 재실행이 완료되지 않았습니다.

## 포함 범위

- 원 형태거리 행렬과 보관 분자거리 행렬
- 저장 상태, feature 방향 교정, 최소교정, 표준화 sequence-set의 분자거리 행렬
- accession 교정표와 확인된 입력 문제
- focal 여섯 종의 `rbcL`, `matK` sequence set
- 선별 alignment, treefile, consensus tree와 IQ-TREE 요약
- 원 결과와 교정 결과를 비교하는 그림
- 저장 행렬의 exact Mantel 수치를 검산하는 공개 준비용 스크립트

## 제외 범위

NCBI 다운로드 cache, 전체 17종 GenBank record, 사용자 경로가 든 실행 JSON·로그, MAFFT·IQ-TREE 바이너리와 Codex가 만든 전체 pipeline은 공개본에서 제외했습니다. 따라서 이 폴더만으로 NCBI 다운로드부터 tree 추론까지를 한 명령으로 재실행할 수는 없습니다.

`code/verify_mantel.py --check`는 tree를 다시 추론하지 않습니다. 저장된 형태·분자거리 행렬에서 720개 exact taxon-label permutation을 다시 계산하고 `data/mantel_all_sensitivity.csv`와 일치하는지만 확인합니다.

