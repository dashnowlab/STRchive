import SchemaForm from "@/components/SchemaForm";
import schema from "~/STRchive-loci.schema.json";

/** add some metadata fields */
schema.properties = {
  "edit-title": {
    type: "string",
    title: "Title",
    description: "Succinct title describing these changes",
    section: "Edit",
    examples: ["Update locus information", "Fix allele details"],
  },
  "edit-description": {
    type: "string",
    title: "Description",
    description: "Detailed summary and justification of these changes",
    multiline: true,
    section: "Edit",
  },
  ...schema.properties,
};

schema.required.push("edit-title", "edit-description");

/** new/edit locus form */
const EditForm = ({ locus }) => (
  <SchemaForm
    schema={schema}
    data={locus}
    sections={[
      "Edit",
      "Overview",
      "IDs",
      "Locus",
      "Alleles",
      "Disease",
      "References",
      "Additional Literature",
    ]}
    onSubmit={console.info}
  />
);

export default EditForm;
