import { useMemo } from "react";
import { FaXmark } from "react-icons/fa6";
import { LuFeather, LuSend } from "react-icons/lu";
import { cloneDeep, isEqual, omitBy, startCase } from "lodash-es";
import { useLocalStorage } from "@reactuses/core";
import { createPR } from "@/api/pr";
import Alert from "@/components/Alert";
import Button from "@/components/Button";
import { contactSchema } from "@/components/ContactForm";
import Form from "@/components/Form";
import Heading from "@/components/Heading";
import Link from "@/components/Link";
import SchemaForm from "@/components/SchemaForm";
import { repo } from "@/layouts/meta";
import { useQuery } from "@/util/hooks";
import { shortenUrl } from "@/util/string";
import loci from "~/STRchive-loci.json";
import schema from "~/STRchive-loci.schema.json";

/** add extra fields for edit metadata */
schema.properties = {
  "edit-name": {
    section: "Edit",
    ...contactSchema.name,
    type: ["string", "null"],
    default: "",
  },
  "edit-username": {
    section: "Edit",
    ...contactSchema.username,
    type: ["string", "null"],
    pattern: "^@",
    default: "",
  },
  "edit-email": {
    section: "Edit",
    ...contactSchema.email,
    type: ["string", "null"],
    default: "",
  },
  "edit-title": {
    section: "Edit",
    title: "Edit Title",
    description: "Succinct title describing these changes",
    examples: ["Fix mechanism details", "Update disease onset information"],
    type: "string",
    default: null,
  },
  "edit-description": {
    section: "Edit",
    title: "Edit Description",
    description:
      "Summary of changes, justification for changes, uncertainty in literature, or anything else we should know for review. Please be detailed. Provide at least 2-3 sentences.",
    examples: [
      "Currently, the disease mechanism details cite a recently retracted paper doi:123456. This edit corrects the reference and updates...",
    ],
    multiline: true,
    type: "string",
    default: null,
  },
  ...schema.properties,
};

schema.required.push("edit-title", "edit-description");

/** new/edit locus form */
const EditForm = ({ heading, locus }) => {
  /** confirm with user before leaving page */
  // window.onbeforeunload = () => "";

  /** unique storage key for this page and form */
  const storageKey = `edit-locus-${locus?.id ?? "new"}`;

  /** form data state */
  let [data, setData] = useLocalStorage(storageKey, {
    ["edit-title"]: null,
    ["edit-description"]: null,
    ...cloneDeep(locus),
  });

  /** was data loaded from storage */
  const storageExists = useMemo(() => {
    const fromStorage = window.localStorage.getItem(storageKey);
    /** if saved draft exists and is different from initial data */
    return fromStorage && !isEqual(JSON.parse(fromStorage), data);
  }, []);

  /** submission query */
  const {
    query: submit,
    data: response,
    status,
  } = useQuery(async () => {
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

    /** remove edit metadata */
    data = omitBy(cloneDeep(data), (value, key) => key.startsWith("edit-"));

    /** merge with locus data */
    data = cloneDeep(loci).map((locus) =>
      locus.id === data.id ? data : locus,
    );

    /** pr files to change */
    const files = [
      {
        path: "data/STRchive-loci.json",
        content: JSON.stringify(data, null, 2),
      },
    ];

    return await createPR({
      owner: "dashnowlab",
      repo: "STRchive",
      branch: locus.id,
      title: data["edit-title"],
      body,
      files,
      labels: ["locus-edit"],
    });
  });

  return (
    <Form onSubmit={submit}>
      <section>
        <Heading level={1}>
          <LuFeather />
          {heading}
        </Heading>

        <Alert type="info">
          Every suggestion is reviewed by our team before inclusion in STRchive.
          Please enter as much accurate information as possible.
        </Alert>

        {storageExists && (
          <div className="row">
            Loaded saved draft
            <Button
              design="plain"
              onClick={() => {
                window.localStorage.removeItem(storageKey);
                window.location.reload();
              }}
            >
              <span>Forget</span>
              <FaXmark />
            </Button>
          </div>
        )}
      </section>

      <SchemaForm
        schema={schema}
        data={data}
        onChange={setData}
        sections={[
          "Edit",
          "Overview",
          "Disease",
          "Locus",
          "Alleles",
          "IDs",
          "References",
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

        <Button type="submit" design="bubble" disabled={status === "success"}>
          <LuSend />
          <span>Submit</span>
        </Button>
      </section>
    </Form>
  );
};

export default EditForm;
