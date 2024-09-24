WITH comments AS (
  SELECT 
    PaperId,
    Doi,
    TRIM(COALESCE(extracted_title1, extracted_title2, extracted_title3)) AS extracted_title
  FROM (
    SELECT 
      PaperId,
      Doi,
      OriginalTitle,
      JournalId,
      EstimatedCitation,
      REGEXP_EXTRACT(OriginalTitle, r'["“`]{1,2}([^"”\'\']+)["”\'\']{1,2}') AS extracted_title1,
      REGEXP_EXTRACT(OriginalTitle, r'Comment on: (.+)$') AS extracted_title2,
      REGEXP_EXTRACT(OriginalTitle, r'Comment on (.+)$') AS extracted_title3
    FROM ccnr-success.mag.Papers
    WHERE JournalId = 77047749
      AND LOWER(OriginalTitle) LIKE "comment%"
      AND Year >= 2000
      AND Year <= 2020
  )
)

SELECT 
  o.Doi as original_doi,
  c.Doi as letter_doi
FROM comments c 
INNER JOIN (
  SELECT 
    PaperId, Doi, OriginalTitle
  FROM ccnr-success.mag.Papers p 
  WHERE JournalId = 77047749
) o 
ON EDIT_DISTANCE(LOWER(c.extracted_title), LOWER(o.OriginalTitle)) < 5
WHERE o.Doi IS NOT NULL 
  AND c.Doi IS NOT NULL
  AND c.Doi != ''
  AND o.Doi != ''
