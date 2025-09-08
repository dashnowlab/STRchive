import SchemaForm from "@/components/SchemaForm";
import schema from "~/STRchive-loci.schema.json";

/** add some metadata fields */
schema.properties = {
  "edit-name": {
    section: "Edit",
    title: "Name",
    placeholder: "Your Name",
    description: "Optional. So we know who you are.",
    type: ["string", "null"],
    default: "",
  },
  "edit-username": {
    section: "Edit",
    title: "GitHub Username",
    placeholder: "@username",
    description: "Optional. So we can tag you.",
    type: ["string", "null"],
    pattern: "^@",
    default: "",
  },
  "edit-email": {
    section: "Edit",
    title: "Email",
    placeholder: "your.name@email.com",
    description: "Optional. So we can contact you directly if needed.",
    type: ["string", "null"],
    default: "",
  },
  "edit-title": {
    section: "Edit",
    title: "Title",
    description: "Succinct title describing these changes",
    examples: ["Fix allele details", "Update locus information"],
    type: "string",
    default: null,
  },
  "edit-description": {
    section: "Edit",
    title: "Description",
    description: "High-level summary and justification of these changes",
    multiline: true,
    type: "string",
    default: null,
  },
  ...schema.properties,
};

schema.required.push("edit-title", "edit-description");

/** new/edit locus form */
const EditForm = ({ locus }) => {
  /** confirm with user before leaving page */
  // window.onbeforeunload = () => "";

  return (
    <SchemaForm
      schema={schema}
      data={{ ["edit-title"]: null, ["edit-description"]: null, ...locus }}
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
};

export default EditForm;
