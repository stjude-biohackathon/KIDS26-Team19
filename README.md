# Team 19 Project

This repository is a starting point for a three-day team project. This repository is populated with a starting template for team organization and planning. Use it to plan, build, and document work. Please adjust this repository to suit the needs of your team.

> **Team leads:** Start with the [team lead checklist](project-management/CHECKLIST.md) before the event or during your first team meeting.

## Project Profile

- **Project name:** Integrating Molecular and Clinical Trial Data to Predict AML Patient Outcomes for Novel Therapies
- **Question, problem, or opportunity:** Harmonizing unstandardized GEO AML data to extract treatment information for prognostic scoring
- **Data, inputs, or evidence:** Gene Expression Omnibus
- **Expected output:** R Shiny Application
- **Tools and stack:** R and Python
- **Team lead:** Nobel Makonnen [SmartOval]
- **Team members and roles:** [Link to `project-management/team.md`]
- **Communication:** Team19 Slack Channel

Naming the tools and stack early helps the team lead create useful roles and divide work realistically. It is fine to revise this section as the project develops.

## Vision and Mission

- **Vision:** To ease discovery of promising therapies for researchers and clinicians
- **Mission:** Collaborate and Incorporate various viewpoints into creating a novel tool for AML researchers

## About

Roughly 30-40% of children with AML given standard therapy will die within 3 years of their diagnosis.  Several prognostic molecular signatures have been developed for pediatric AML.  My team developed a prognostic molecular signature for pediatric AML, found that decitabine tretment greatly improved this marker, and thereby estimated that adding decitabine to standard AML therapy could cut the poor outcome rate in half.  This provided the scientific rationale for the AML16 clinical trial which confirmed this prediction.  In this project, we will mine publicly available datasets to estimate drugs’ impacts on prognostic markers to identify additional drugs that hold promise to improve outcomes for pediatric AML.  

## Roadmap and Milestones

| When | Focus | Expected outcome |
| --- | --- | --- |
| Day 1 | Introductions, Assigning roles, and beginning to build functions to query and download studies from GEO |
| Day 2 | Building a database to house the studies and working on an R shiny app to house the tool |
| Day 3 | Finishing touches and validating that the tool works |

The goal is not a perfect production system. The goal is a clear, honest, useful result that the team can explain and others can build on.

## Project Setup
_This section contains draft content. Please check back towards the end of the project._

To install python dependencies, run from the root:
```
pip install -r parsing/requirement.txt
```


