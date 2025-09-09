import { useRef } from "react";
import { LuSend } from "react-icons/lu";
import { startCase } from "lodash-es";
import { createPR } from "@/api/pr";
import Alert from "@/components/Alert";
import Button from "@/components/Button";
import Form from "@/components/Form";
import Link from "@/components/Link";
import SchemaForm from "@/components/SchemaForm";
import { repo } from "@/layouts/meta";
import { useQuery } from "@/util/hooks";
import { shortenUrl } from "@/util/string";
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
  const ref = useRef();

  /** confirm with user before leaving page */
  // window.onbeforeunload = () => "";

  /** submission query */
  const {
    query: submit,
    data: response,
    status,
  } = useQuery(async () => {
    /** get current form data */
    const { data } = ref.current;

    /** pr branch name */
    const branch = locus.id;

    /** pr title */
    const title = data["edit-title"];

    /** pr body */
    const name = data["edit-name"];
    const username = data["edit-username"];
    const email = data["edit-email"];
    const description = data["edit-description"];
    const body = [{ name, username, email }, { description }]
      .map((group) =>
        Object.entries(group)
          .map(([key, value]) => [
            `**${startCase(key)}**`,
            value?.trim() ? value.trim() : "\\-",
          ])
          .flat()
          .join("\n"),
      )
      .join("\n\n");

    /** pr files to change */
    const files = [
      {
        path: "data/STRchive-loci.json",
        content: JSON.stringify(data, null, 2),
      },
    ];

    /** pr labels */
    const labels = ["locus-edit"];

    return await createPR({ branch, title, body, files, labels });
  });

  return (
    <Form onSubmit={submit}>
      {/* maintain correct section coloring */}
      <span></span>

      <SchemaForm
        ref={ref}
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
        onSubmit={(data) => submit(data)}
      />

      <section>
        <Alert type={status || "info"}>
          {status === "" && (
            <>
              This will make a <strong>public</strong> pull request on{" "}
              <Link to={repo}>our GitHub</Link>. You'll get a link once it's
              created.
            </>
          )}
          {startCase(status)}{" "}
          {response && (
            <Link to={response.link}>{shortenUrl(response.link)}</Link>
          )}
        </Alert>

        <Button type="submit" design="bubble">
          <LuSend />
          <span>Submit</span>
        </Button>
      </section>
    </Form>
  );
};

export default EditForm;
